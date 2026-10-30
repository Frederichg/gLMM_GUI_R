"""
R-backed GLMM analyzer GUI (first robust version).

This app keeps the familiar Tkinter workflow but delegates model fitting,
diagnostics, emmeans, and report generation to R via Rscript.


24 avril 2026, add save as the input file name instead of timestamp to make it easier to track outputs from different datasets.
9 juin 2026, add to output all the Emmeans for each VD by VI.
12 juin 2026, change default choices and add percent scales on 1 or 100
22 juin 2026, add nesting and continuous to the IV. Display the R formula.
02 juillet 2026, added editable R equation
6 juillet 2026, added R equation in the output. Added Timestamp and Computer ID to the output. Many corrections to the reports. 
8 juillet 2026, gaussian LMM pour les données continues non-skew. 
                cbind for the minomial data.
31 juillet 2026, correct some bugs with the formula
30 oct 2026, recreate the 'easy' published version
"""

from __future__ import annotations

import json
import re
import socket
import subprocess
import threading
from datetime import datetime
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import pandas as pd


APP_DIR = Path(__file__).resolve().parent
CONFIG_PATH = APP_DIR / ".r_glmm_config.json"
RUN_ANALYSIS_SCRIPT = APP_DIR / "run_analysis.R"
PROGRAM_NAME = Path(__file__).resolve().name

IV_ROLE_VALUES = ["fixed_effect", "random_intercept", "random_effect", "random_slope"]
IV_SCOPE_VALUES = ["within_subject", "between_subject"]
IV_STRUCTURE_VALUES = ["nothing", "crossed", "nested", "continuous"]
# Visible DV type choices (user-facing). We map these to internal type keys when collecting payload.
DV_TYPE_VALUES = [
    "Poisson (count)",
    "Gamma (reaction time)",
    "Gaussian (continuous non-skewed)",
    "Binomial (proportion)",
    "Beta (100% = 1)",
    "Percent 100% = 100",
]

# Map user-visible DV type labels (lower-cased) to internal type keys or special handling markers.
DISPLAY_TO_INTERNAL_DV = {
    "poisson (count)": "counts",
    "gamma (reaction time)": "rt",
    "gaussian (continuous non-skewed)": "gaussian",
    "binomial (proportion)": "binomial_cbind",
    "beta (100% = 1)": "proportion",
    "counts": "counts",
    "reaction time": "rt",
    "continuous non-skewed": "gaussian",
    "proportion (binomial - success/failure counts)": "binomial_cbind",
    "proportion_without_non_visible": "proportion_without_non_visible",
    "proportion_multinomial": "proportion_multinomial",
    "percent 100% = 1": "proportion",
    "percent 100% = 100": "proportion_percent100",
}


class RGLMMAnalyzerGUI:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("GLMM Analyzer (R Backend)")
        self.root.geometry("1000x1000")

        self.data_path: Path | None = None
        self.data_df: pd.DataFrame | None = None
        self.sheet_name: str | None = None
        self.last_output_dir: Path | None = None
        self._column_list: list[str] = []  # Store full column list for cbind pair combos

        self._build_gui()

    def _build_gui(self) -> None:
        main = ttk.Frame(self.root)
        main.pack(fill="both", expand=True)

        top = ttk.Frame(main)
        top.pack(fill="x", padx=10, pady=10)

        self._build_file_section(top)
        self._build_structure_section(top)
        self._build_actions(top)

        bottom = ttk.LabelFrame(main, text="Execution Log", padding=10)
        bottom.pack(fill="both", expand=True, padx=10, pady=(0, 10))

        self.log_text = tk.Text(bottom, height=22, wrap="word")
        self.log_text.pack(fill="both", expand=True)

    def _build_file_section(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Data File", padding=10)
        frame.pack(fill="x", pady=5)

        self.file_label = ttk.Label(frame, text="No file selected", foreground="gray")
        self.file_label.pack(side="left", fill="x", expand=True)

        ttk.Button(frame, text="Browse Excel File...", command=self._browse_file).pack(side="right")

    def _build_structure_section(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Analysis Structure", padding=10)
        frame.pack(fill="x", pady=5)

        # Sheet selector
        sheet_row = ttk.Frame(frame)
        sheet_row.pack(fill="x", pady=4)
        ttk.Label(sheet_row, text="Sheet:", width=16).pack(side="left")
        self.sheet_combo = ttk.Combobox(sheet_row, state="readonly", width=40)
        self.sheet_combo.pack(side="left", padx=5)
        self.sheet_combo.bind("<<ComboboxSelected>>", lambda _: self._load_columns_for_selected_sheet())

        # Animal ID
        id_row = ttk.Frame(frame)
        id_row.pack(fill="x", pady=4)
        ttk.Label(id_row, text="Animal ID:", width=16).pack(side="left")
        self.id_combo = ttk.Combobox(id_row, state="readonly", width=40)
        self.id_combo.pack(side="left", padx=5)

        # VIs + role/scope
        vi_box = ttk.LabelFrame(frame, text="Independent Variables (up to 4)", padding=6)
        vi_box.pack(fill="x", pady=6)
        help_row = ttk.Frame(vi_box)
        help_row.pack(fill="x", pady=(0, 4))
        ttk.Label(
            help_row,
            text="Choisissez pour chaque VI son role dans le modele et sa structure sujet.",
            foreground="gray",
            font=("TkDefaultFont", 8),
        ).pack(side="left")
        ttk.Button(help_row, text="?", width=3, command=self._show_vi_help).pack(side="left", padx=6)

        self.vi_combos: list[ttk.Combobox] = []
        self.vi_role_combos: list[ttk.Combobox] = []
        self.vi_scope_combos: list[ttk.Combobox] = []
        self.vi1_structure_combo: ttk.Combobox | None = None
        self.vi1_cross_combo: ttk.Combobox | None = None
        self.vi1_cross_label: ttk.Label | None = None
        self.vi1_nested_combos: list[ttk.Combobox] = []
        self.vi1_structure_info: ttk.Label | None = None

        for i in range(4):
            row = ttk.Frame(vi_box)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"VI {i + 1}:", width=8).pack(side="left")
            vi_combo = ttk.Combobox(row, state="readonly", width=28)
            vi_combo.pack(side="left", padx=4)
            vi_combo.bind("<<ComboboxSelected>>", lambda _e: self._update_formula_text())
            self.vi_combos.append(vi_combo)

            ttk.Label(row, text="Role:", width=6).pack(side="left", padx=(10, 2))
            vi_role = ttk.Combobox(
                row,
                state="readonly",
                width=16,
                values=IV_ROLE_VALUES,
            )
            vi_role.set("fixed_effect")
            vi_role.pack(side="left")
            vi_role.bind("<<ComboboxSelected>>", lambda _e: self._update_formula_text())
            self.vi_role_combos.append(vi_role)

            ttk.Label(row, text="Scope:", width=7).pack(side="left", padx=(10, 2))
            vi_scope = ttk.Combobox(
                row,
                state="readonly",
                width=16,
                values=IV_SCOPE_VALUES,
            )
            # Default to between-subject as requested
            vi_scope.set("between_subject")
            vi_scope.pack(side="left")
            vi_scope.bind("<<ComboboxSelected>>", lambda _e: self._update_formula_text())
            self.vi_scope_combos.append(vi_scope)

            if i == 0:
                self.vi1_structure_combo = ttk.Combobox(
                    row,
                    state="readonly",
                    width=14,
                    values=IV_STRUCTURE_VALUES,
                )
                self.vi1_structure_combo.set("nothing")
                self.vi1_structure_combo.bind("<<ComboboxSelected>>", lambda _e: self._update_vi1_structure_state())

        self.vi1_structure_frame = ttk.Frame(vi_box)

        self.vi1_cross_label = ttk.Label(self.vi1_structure_frame, text="Cross with:", width=16)
        self.vi1_cross_combo = ttk.Combobox(self.vi1_structure_frame, state="readonly", width=28)
        self.vi1_cross_combo.bind("<<ComboboxSelected>>", lambda _e: self._update_formula_text())

        for _ in range(3):
            nested_combo = ttk.Combobox(self.vi1_structure_frame, state="readonly", width=20)
            nested_combo.bind("<<ComboboxSelected>>", lambda _e: self._update_formula_text())
            self.vi1_nested_combos.append(nested_combo)

        self.vi1_structure_info = ttk.Label(self.vi1_structure_frame, text="", foreground="orange")

        formula_frame = ttk.LabelFrame(frame, text="R Formula (Edit for Advanced Use)", padding=6)
        formula_frame.pack(fill="both", expand=False, pady=5)

        formula_controls = ttk.Frame(formula_frame)
        formula_controls.pack(fill="x", pady=(0, 4))
        ttk.Button(formula_controls, text="Apply Formula from GUI", command=self._apply_formula_from_gui).pack(side="left", padx=2)
        ttk.Button(formula_controls, text="Reset to Default", command=self._reset_formula).pack(side="left", padx=2)
        ttk.Label(formula_controls, text="(Single-line formula: edit to add interactions like '*' or modify)", foreground="gray", font=("TkDefaultFont", 8)).pack(side="left", padx=8)

        formula_warning = ttk.Frame(formula_frame)
        formula_warning.pack(fill="x", pady=(0, 4))
        ttk.Label(formula_warning, text="⚠ Column names in formula must match data EXACTLY (case-sensitive, e.g. 'phenotype' not 'Phenotype')", foreground="darkred", font=("TkDefaultFont", 8)).pack(side="left", padx=2)

        self.formula_text = tk.Text(formula_frame, height=1, width=100, wrap="none")
        self.formula_text.pack(fill="both", expand=False)
        self.formula_text.bind("<KeyRelease>", lambda _: self._mark_formula_edited())
        self.formula_text.bind("<Return>", lambda e: "break")

        self.formula_var = tk.StringVar(value="DV ~ ... + (1 | Animal_ID)")
        self.formula_edited = False
        self.auto_generated_formula = "DV ~ ... + (1 | Animal_ID)"
        self.formula_text.insert("1.0", self.formula_var.get())

        self.vd_combos: list[ttk.Combobox] = []
        self.vd_type_combos: list[ttk.Combobox] = []
        self.vd_cbind_pair_combos: list[ttk.Combobox] = []
        self.vd_cbind_pair_labels: list[ttk.Label] = []

        self._update_vi1_structure_state()

        # VDs + type (+ non-visible for dedicated proportion mode)
        vd_box = ttk.LabelFrame(frame, text="Dependent Variables (up to 10)", padding=6)
        vd_box.pack(fill="x", pady=6)

        cycle_row = ttk.Frame(vd_box)
        ttk.Label(cycle_row, text="Observation cycle per row:", width=26).pack(side="left")
        self.obs_cycle_var = tk.StringVar(value="30")
        self.obs_cycle_spin = ttk.Spinbox(cycle_row, from_=1, to=10000, textvariable=self.obs_cycle_var, width=8)
        self.obs_cycle_spin.pack(side="left")
        ttk.Label(
            cycle_row,
            text="Used by 'proportion_without_non_visible' (default 30).",
            foreground="gray",
            font=("TkDefaultFont", 8),
        ).pack(side="left", padx=8)
        non_visible_row = ttk.Frame(vd_box)
        ttk.Label(non_visible_row, text="Non-visible column:", width=26).pack(side="left")
        self.global_nonvisible_combo = ttk.Combobox(non_visible_row, state="disabled", width=28)
        self.global_nonvisible_combo.pack(side="left")
        ttk.Label(
            non_visible_row,
            text="Shared by all DVs using proportion_without_non_visible / proportion_multinomial.",
            foreground="gray",
            font=("TkDefaultFont", 8),
        ).pack(side="left", padx=8)

        for i in range(10):
            row = ttk.Frame(vd_box)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"VD {i + 1}:", width=8).pack(side="left")
            vd_combo = ttk.Combobox(row, state="readonly", width=24)
            vd_combo.pack(side="left", padx=4)
            vd_combo.bind("<<ComboboxSelected>>", lambda _e: self._update_cbind_pair_visibility())
            self.vd_combos.append(vd_combo)

            ttk.Label(row, text="Type:", width=6).pack(side="left", padx=(10, 2))
            vd_type = ttk.Combobox(
                row,
                state="readonly",
                width=30,
                values=DV_TYPE_VALUES,
            )
            vd_type.set("Poisson (count)")
            vd_type.pack(side="left")
            self.vd_type_combos.append(vd_type)

            cbind_pair_label = ttk.Label(row, text="Pair with:", width=10)
            self.vd_cbind_pair_labels.append(cbind_pair_label)

            cbind_pair_combo = ttk.Combobox(row, state="readonly", width=20)
            self.vd_cbind_pair_combos.append(cbind_pair_combo)

            vd_type.bind("<<ComboboxSelected>>", lambda _e, idx=i: self._on_vd_type_changed(idx))

    def _build_actions(self, parent: ttk.Frame) -> None:
        frame = ttk.LabelFrame(parent, text="Actions", padding=10)
        frame.pack(fill="x", pady=5)

        self.status_var = tk.StringVar(value="Idle")
        ttk.Label(frame, textvariable=self.status_var, foreground="blue").pack(side="left")

        self.analyze_button = ttk.Button(frame, text="Analyze", command=self._start_analysis)
        self.analyze_button.pack(side="right", padx=5)

        ttk.Button(frame, text="Open Last Output", command=self._open_last_output).pack(side="right", padx=5)
        ttk.Button(frame, text="Set Rscript Path", command=self._choose_rscript_path).pack(side="right", padx=5)

        # Label showing the currently saved Rscript path
        saved_rscript = self._load_config().get("rscript_path", "").strip()
        display = saved_rscript if saved_rscript else "Rscript (auto-detect from PATH)"
        self.rscript_label = ttk.Label(frame, text=f"Rscript: {display}", foreground="gray", font=("TkDefaultFont", 8))
        self.rscript_label.pack(side="right", padx=(0, 10))

    def _apply_formula_from_gui(self) -> None:
        """Apply the auto-generated formula from GUI selections to the text widget."""
        self._update_formula_text()
        self.formula_text.delete("1.0", "end")
        self.formula_text.insert("1.0", self.formula_var.get())
        self.formula_edited = False
        self._log("Formula updated from GUI selections.")

    def _reset_formula(self) -> None:
        """Reset the formula to the default auto-generated formula."""
        self.formula_edited = False
        self._update_formula_text()
        self.formula_text.delete("1.0", "end")
        self.formula_text.insert("1.0", self.formula_var.get())
        self.auto_generated_formula = self.formula_var.get()
        self._log("Formula reset to default.")

    def _mark_formula_edited(self) -> None:
        """Mark that the formula has been manually edited by the user."""
        self.formula_edited = True

    def _browse_file(self) -> None:
        path = filedialog.askopenfilename(
            title="Select Excel file",
            filetypes=[("Excel", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if not path:
            return

        try:
            self.data_path = Path(path)
            xl = pd.ExcelFile(self.data_path)
            sheets = xl.sheet_names
            self.sheet_combo["values"] = sheets
            if sheets:
                self.sheet_combo.set(sheets[0])
                self._load_columns_for_selected_sheet()

            self.file_label.config(text=str(self.data_path), foreground="black")
            self._log(f"Loaded workbook: {self.data_path}")
        except Exception as exc:
            messagebox.showerror("Error", f"Unable to read workbook:\n{exc}")

    def _clear_column_selections(self) -> None:
        self.id_combo.set("")
        for combo in self.vi_combos:
            combo.set("")
        for combo in self.vd_combos:
            combo.set("")
        self.global_nonvisible_combo.set("")
        if self.vi1_cross_combo is not None:
            self.vi1_cross_combo.set("")
        for nested_combo in self.vi1_nested_combos:
            nested_combo.set("")
        self._update_global_nonvisible_state()
        self._update_formula_text()

    def _load_columns_for_selected_sheet(self) -> None:
        if not self.data_path:
            return

        sheet = self.sheet_combo.get().strip()
        if not sheet:
            return

        try:
            self.sheet_name = sheet
            self.data_df = pd.read_excel(self.data_path, sheet_name=sheet, engine="openpyxl")
            self.data_df.columns = [str(c).strip() for c in self.data_df.columns]
            columns = [""] + [str(c) for c in self.data_df.columns]
            self._column_list = columns

            self._clear_column_selections()
            self.id_combo["values"] = columns
            for combo in self.vi_combos:
                combo["values"] = columns
            for combo in self.vd_combos:
                combo["values"] = columns
            for combo in self.vd_cbind_pair_combos:
                combo["values"] = columns
            self.global_nonvisible_combo["values"] = columns
            if self.vi1_cross_combo is not None:
                self.vi1_cross_combo["values"] = columns
            for nested_combo in self.vi1_nested_combos:
                nested_combo["values"] = columns
            self._update_global_nonvisible_state()
            self._update_cbind_pair_visibility()

            self._log(f"Selected sheet: {sheet} ({len(self.data_df)} rows, {len(self.data_df.columns)} columns)")
        except Exception as exc:
            messagebox.showerror("Error", f"Unable to read selected sheet:\n{exc}")

    def _show_vi_help(self) -> None:
        help_lines = [
            "Guide rapide (francais):",
            "",
            "1) fixed_effect:",
            "Si vous voulez tester si une variable influence votre comportement, choisissez fixed_effect.",
            "Exemple: 'Sexe' (M/F) -> between_subject + fixed_effect.",
            "",
            "2) random_intercept:",
            "Utilisez random_intercept pour une source de variabilite de groupe comme le numero de cage, le lot ou la cohorte.",
            "Exemple: numero de cage -> between_subject + random_intercept.",
            "Cela ajoute un intercept aleatoire du type (1|cage).",
            "",
            "3) random_effect:",
            "random_effect peut etre utilise comme categorie generale si vous savez qu'une variable doit etre modelee en aleatoire mais sans precision supplementaire.",
            "Dans la pratique, pour une cage ou une cohorte, privilegiez random_intercept.",
            "Exemple: 'Cohorte' -> between_subject + random_effect.",
            "",
            "4) random_slope:",
            "Utilisez random_slope pour une variable qui varie dans le sujet et dont la pente peut changer selon l'animal.",
            "Exemple: traitement repete dans le temps -> within_subject + random_slope.",
            "",
            "5) within_subject vs between_subject:",
            "within_subject: la variable change a l'interieur du meme animal.",
            "between_subject: la variable differencie les animaux entre eux.",
            "Exemple within_subject: AM/PM (mesure repetee sur le meme animal).",
            "Exemple between_subject: Sexe ou Cage (caracteristique de l'animal/groupe).",
        ]
        messagebox.showinfo("Aide VI", "\n".join(help_lines))

    def _update_global_nonvisible_state(self) -> None:
        requires_non_visible = any(
            DISPLAY_TO_INTERNAL_DV.get(combo.get().strip().lower(), combo.get().strip().lower())
            in ("proportion_without_non_visible", "proportion_multinomial")
            for combo in self.vd_type_combos
        )
        if requires_non_visible:
            self.global_nonvisible_combo.configure(state="readonly")
        else:
            self.global_nonvisible_combo.set("")
            self.global_nonvisible_combo.configure(state="disabled")

    def _on_vd_type_changed(self, vd_index: int = 0) -> None:
        self._update_global_nonvisible_state()
        self._update_cbind_pair_visibility()
        self._update_formula_text()

    def _update_cbind_pair_visibility(self) -> None:
        col_values = self._column_list if self._column_list else [""]
        for vd_type_combo, cbind_pair_label, cbind_pair_combo in zip(
            self.vd_type_combos,
            self.vd_cbind_pair_labels,
            self.vd_cbind_pair_combos,
        ):
            dv_display = vd_type_combo.get().strip()
            dv_internal = DISPLAY_TO_INTERNAL_DV.get(dv_display.lower(), dv_display.lower())

            if dv_internal == "binomial_cbind":
                cbind_pair_combo["values"] = col_values
                if not cbind_pair_label.winfo_manager():
                    cbind_pair_label.pack(side="left", padx=(15, 2))
                if not cbind_pair_combo.winfo_manager():
                    cbind_pair_combo.pack(side="left", padx=2)
            else:
                cbind_pair_label.pack_forget()
                cbind_pair_combo.pack_forget()
                cbind_pair_combo.set("")

    def _update_vi1_structure_state(self) -> None:
        structure = self.vi1_structure_combo.get().strip().lower() if self.vi1_structure_combo else "nothing"
        if structure == "crossed":
            self.vi1_cross_label.configure(foreground="black")
            self.vi1_cross_combo.configure(state="readonly")
            for nested_combo in self.vi1_nested_combos:
                nested_combo.configure(state="disabled")
                nested_combo.set("")
            self.vi1_structure_info.configure(text="Use this to add an interaction term with VI1.")
        elif structure == "nested":
            self.vi1_cross_label.configure(foreground="gray")
            self.vi1_cross_combo.configure(state="disabled")
            self.vi1_cross_combo.set("")
            for nested_combo in self.vi1_nested_combos:
                nested_combo.configure(state="readonly")
            self.vi1_structure_info.configure(text="Select up to 3 nesting levels for VI1.")
        elif structure == "continuous":
            self.vi1_cross_label.configure(foreground="gray")
            self.vi1_cross_combo.configure(state="disabled")
            self.vi1_cross_combo.set("")
            for nested_combo in self.vi1_nested_combos:
                nested_combo.configure(state="disabled")
                nested_combo.set("")
            self.vi1_structure_info.configure(text="Continuous IV will be z-transformed; avoid high correlation with other IVs.")
        else:
            self.vi1_cross_label.configure(foreground="gray")
            self.vi1_cross_combo.configure(state="disabled")
            self.vi1_cross_combo.set("")
            for nested_combo in self.vi1_nested_combos:
                nested_combo.configure(state="disabled")
                nested_combo.set("")
            self.vi1_structure_info.configure(text="")
        self._update_formula_text()

    def _update_formula_text(self) -> None:
        dv_parts = []
        processed_cbind = set()
        for i, (vd_combo, vd_type_combo) in enumerate(zip(self.vd_combos, self.vd_type_combos)):
            vd_name = vd_combo.get().strip()
            if not vd_name:
                continue
            if i in processed_cbind:
                continue

            dv_display = vd_type_combo.get().strip()
            dv_internal = DISPLAY_TO_INTERNAL_DV.get(dv_display.lower(), dv_display.lower())

            if dv_internal == "binomial_cbind":
                pair_name = self.vd_cbind_pair_combos[i].get().strip()
                if pair_name:
                    dv_parts.append(f"cbind(`{vd_name}`, `{pair_name}`)")
                    for j, other_combo in enumerate(self.vd_combos):
                        if other_combo.get().strip() == pair_name:
                            processed_cbind.add(j)
                else:
                    dv_parts.append(f"`{vd_name}`")
            else:
                dv_parts.append(f"`{vd_name}`")

        dv_part = "DV"

        iv_terms = []
        for i, vi_combo in enumerate(self.vi_combos):
            name = vi_combo.get().strip()
            if not name:
                continue
            if i == 0 and self.vi1_structure_combo is not None:
                structure = self.vi1_structure_combo.get().strip().lower()
                if structure == "crossed":
                    cross_with = self.vi1_cross_combo.get().strip() if self.vi1_cross_combo else ""
                    if cross_with:
                        iv_terms.append(f"`{name}` * `{cross_with}`")
                    else:
                        iv_terms.append(f"`{name}`")
                    continue
                if structure == "nested":
                    nested = [combo.get().strip() for combo in self.vi1_nested_combos if combo.get().strip()]
                    if nested:
                        nested_terms = " %in% ".join([f"`{n}`" for n in nested])
                        iv_terms.append(f"`{name}` %in% {nested_terms}")
                    else:
                        iv_terms.append(f"`{name}`")
                    continue
                if structure == "continuous":
                    iv_terms.append(f"scale(`{name}`)")
                    continue
            iv_terms.append(f"`{name}`")

        formula = " + ".join(iv_terms) if iv_terms else "..."
        animal_id_col = self.id_combo.get().strip()
        animal_id_expr = f"(1 | `{animal_id_col}`)" if animal_id_col else "(1 | Animal_ID)"
        new_formula = f"DV ~ {formula} + {animal_id_expr}"
        self.formula_var.set(new_formula)
        self.auto_generated_formula = new_formula
        if not self.formula_edited:
            self.formula_text.delete("1.0", "end")
            self.formula_text.insert("1.0", new_formula)

    def _choose_rscript_path(self) -> None:
        cfg = self._load_config()
        current = cfg.get("rscript_path", "").strip()
        initial_dir = str(Path(current).parent) if current and Path(current).exists() else "/"

        path = filedialog.askopenfilename(
            title="Select Rscript executable",
            initialdir=initial_dir,
            filetypes=[("Rscript", "Rscript.exe"), ("Executable", "*.exe"), ("All", "*.*")],
        )
        if not path:
            return

        cfg["rscript_path"] = path
        self._save_config(cfg)
        self.rscript_label.config(text=f"Rscript: {path}", foreground="black")
        self._log(f"Saved Rscript path: {path}")

    def _start_analysis(self) -> None:
        if not self.data_path or self.data_df is None or not self.sheet_name:
            messagebox.showerror("Missing data", "Please select an Excel file and sheet first.")
            return

        if not self.id_combo.get().strip():
            messagebox.showerror("Missing Animal ID", "Please select an Animal ID column.")
            return

        self._update_formula_text()
        output_dir = APP_DIR / "outputs" / self.data_path.stem
        output_dir.mkdir(parents=True, exist_ok=True)
        self.last_output_dir = output_dir

        try:
            payload = self._collect_payload()
        except ValueError as exc:
            messagebox.showerror("Missing input", str(exc))
            return

        config_path = output_dir / "analysis_config.json"
        config_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

        self._set_busy(True)
        self._log("Starting R analysis...")

        threading.Thread(
            target=self._run_r_subprocess,
            args=(config_path, output_dir),
            daemon=True,
        ).start()

    def _normalize_column_name(self, column_name: str) -> str:
        if column_name is None:
            return ""
        return re.sub(r"\s+", " ", str(column_name).strip()).casefold()

    def _resolve_column_name(self, column_name: str) -> str:
        if not column_name:
            return ""

        exact_name = str(column_name).strip()
        if exact_name in self.data_df.columns:
            return exact_name

        normalized_target = self._normalize_column_name(exact_name)
        for actual_name in self.data_df.columns:
            if self._normalize_column_name(actual_name) == normalized_target:
                return str(actual_name)
        return ""

    def _collect_payload(self) -> dict:
        if self.data_df is None:
            raise ValueError("No sheet data is currently loaded.")

        available_columns = {str(col) for col in self.data_df.columns}
        id_col = self.id_combo.get().strip()

        def _ensure_column_exists(column_name: str, label: str) -> None:
            if column_name:
                resolved_name = self._resolve_column_name(column_name)
                if not resolved_name:
                    raise ValueError(
                        f"Selected {label} column '{column_name}' is not present in sheet '{self.sheet_name}'."
                        f" Please choose a column from the current sheet."
                    )
                return

        try:
            observation_cycle = int(self.obs_cycle_var.get().strip())
        except ValueError as exc:
            raise ValueError("Observation cycle must be an integer (e.g., 30).") from exc

        if observation_cycle <= 0:
            raise ValueError("Observation cycle must be > 0.")

        resolved_id_col = self._resolve_column_name(id_col) if id_col else ""
        if not resolved_id_col:
            raise ValueError(f"Selected Animal ID column '{id_col}' is not present in sheet '{self.sheet_name}'.")
        id_col = resolved_id_col
        _ensure_column_exists(id_col, "Animal ID")

        vi_items: list[dict] = []
        for i, (vi_combo, vi_role_combo, vi_scope_combo) in enumerate(zip(
            self.vi_combos,
            self.vi_role_combos,
            self.vi_scope_combos,
        )):
            vi = vi_combo.get().strip()
            if vi:
                resolved_vi = self._resolve_column_name(vi)
                if not resolved_vi:
                    raise ValueError(f"Selected independent variable column '{vi}' is not present in sheet '{self.sheet_name}'.")
                _ensure_column_exists(resolved_vi, "independent variable")
                scope = vi_scope_combo.get().strip().lower() or "within_subject"
                item = {
                    "name": resolved_vi,
                    "model_role": vi_role_combo.get().strip().lower() or "fixed_effect",
                    "subject_scope": scope,
                    # Backward-compatible field consumed by older R code paths.
                    "type": "within" if scope == "within_subject" else "between",
                }
                if i == 0:
                    structure = self.vi1_structure_combo.get().strip().lower() if self.vi1_structure_combo else "nothing"
                    item["structure"] = structure
                    if structure == "crossed":
                        cross_with = self.vi1_cross_combo.get().strip() if self.vi1_cross_combo else ""
                        resolved_cross_with = self._resolve_column_name(cross_with)
                        if cross_with and not resolved_cross_with:
                            raise ValueError(f"Selected crossed variable '{cross_with}' is not present in sheet '{self.sheet_name}'.")
                        _ensure_column_exists(resolved_cross_with, "crossed variable")
                        item["crossed_with"] = resolved_cross_with
                    elif structure == "nested":
                        nested_cols = [combo.get().strip() for combo in self.vi1_nested_combos if combo.get().strip()]
                        resolved_nested_cols = []
                        for nested_col in nested_cols:
                            resolved_nested_col = self._resolve_column_name(nested_col)
                            if nested_col and not resolved_nested_col:
                                raise ValueError(f"Selected nested variable '{nested_col}' is not present in sheet '{self.sheet_name}'.")
                            resolved_nested_cols.append(resolved_nested_col)
                        item["nested_with"] = resolved_nested_cols
                    elif structure == "continuous":
                        item["continuous"] = True
                vi_items.append(item)

        if not vi_items:
            raise ValueError("Please select at least one independent variable.")

        vd_items: list[dict] = []
        force_multinomial = False
        shared_non_visible_col = self.global_nonvisible_combo.get().strip()
        conversion_made = False
        continuous_iv_z = None
        if vi_items and vi_items[0].get("structure") == "continuous":
            continuous_iv_z = vi_items[0]["name"]
            if self.data_df is None:
                raise ValueError("Internal error: data not loaded for continuous IV conversion.")
            try:
                numeric_series = pd.to_numeric(self.data_df[continuous_iv_z], errors="coerce")
                mean_value = numeric_series.mean()
                std_value = numeric_series.std(ddof=0)
                if std_value == 0 or pd.isna(std_value):
                    raise ValueError("standard deviation is zero or undefined")
                self.data_df[continuous_iv_z] = (numeric_series - mean_value) / std_value
                conversion_made = True
                self._log(f"Z-transformed continuous IV: {continuous_iv_z}")
            except Exception as exc:
                raise ValueError(f"Unable to z-transform continuous IV {continuous_iv_z}: {exc}") from exc

        for vd_combo, vd_type_combo in zip(
            self.vd_combos,
            self.vd_type_combos,
        ):
            vd = vd_combo.get().strip()
            if vd:
                resolved_vd = self._resolve_column_name(vd)
                if not resolved_vd:
                    raise ValueError(f"Selected dependent variable column '{vd}' is not present in sheet '{self.sheet_name}'.")
                _ensure_column_exists(resolved_vd, "dependent variable")
                dv_display = vd_type_combo.get().strip()
                dv_internal = DISPLAY_TO_INTERNAL_DV.get(dv_display.lower(), dv_display.lower())
                if dv_internal == "proportion_multinomial":
                    # Backward-compatible: reuse existing proportion-without-non-visible pipeline.
                    # Multinomial section is toggled globally via run_multinomial.
                    force_multinomial = True
                    dv_internal = "proportion_without_non_visible"

                # If the user selected percent scale where 100% is represented as 100,
                # convert the data column to 0-1 for modeling (divide by 100).
                if dv_internal == "proportion_percent100":
                    if self.data_df is None:
                        raise ValueError("Internal error: data not loaded for conversion.")
                    try:
                        self.data_df[resolved_vd] = self.data_df[resolved_vd].astype(float) / 100.0
                        conversion_made = True
                    except Exception as exc:
                        raise ValueError(f"Unable to convert column {resolved_vd} to 0-1 scale: {exc}") from exc
                    dv_internal = "proportion"

                item: dict[str, str] = {"name": resolved_vd, "type": dv_internal}
                if dv_internal == "proportion_without_non_visible":
                    if not shared_non_visible_col:
                        raise ValueError(
                            "Please select the shared Non-visible column for DVs using proportion_without_non_visible/proportion_multinomial."
                        )
                    resolved_shared_non_visible_col = self._resolve_column_name(shared_non_visible_col)
                    if not resolved_shared_non_visible_col:
                        raise ValueError(f"Selected non-visible column '{shared_non_visible_col}' is not present in sheet '{self.sheet_name}'.")
                    _ensure_column_exists(resolved_shared_non_visible_col, "non-visible column")
                    item["non_visible_column"] = resolved_shared_non_visible_col

                vd_items.append(item)

        if not vd_items:
            raise ValueError("Please select at least one dependent variable.")

        binomial_cbind_indices = [i for i, item in enumerate(vd_items) if item["type"] == "binomial_cbind"]
        if binomial_cbind_indices:
            for idx in binomial_cbind_indices:
                pair_name = self.vd_cbind_pair_combos[idx].get().strip()
                if not pair_name:
                    raise ValueError(f"Dependent variable '{vd_items[idx]['name']}' is marked as Proportion (binomial), but no pair was selected. Please select the failure/complement counts column in the 'Pair with' dropdown.")
                resolved_pair_name = self._resolve_column_name(pair_name)
                if not resolved_pair_name:
                    raise ValueError(f"Selected pair column '{pair_name}' for cbind is not present in sheet '{self.sheet_name}'.")
                pair_idx = None
                for j, item in enumerate(vd_items):
                    if item["name"] == resolved_pair_name:
                        pair_idx = j
                        break
                vd_items[idx]["cbind_pair_with"] = resolved_pair_name
                if pair_idx is not None:
                    vd_items[idx]["cbind_pair_index"] = pair_idx
        input_path = str(self.data_path)
        # If we converted any columns to 0-1 internally, write a temporary Excel copy
        # and point the R analysis to that file instead of the original.
        if conversion_made:
            if not self.last_output_dir:
                raise ValueError("Internal error: output directory unavailable for converted data.")
            converted_path = Path(self.last_output_dir) / f"{self.data_path.stem}_converted.xlsx"
            try:
                self.data_df.to_excel(converted_path, sheet_name=self.sheet_name or "Sheet1", index=False)
                input_path = str(converted_path)
                self._log(f"Wrote converted data to {converted_path}")
            except Exception as exc:
                raise ValueError(f"Unable to write converted data file: {exc}") from exc

        run_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        return {
            "input": input_path,
            "sheet": self.sheet_name,
            "animal_id": id_col,
            "iv": vi_items,
            "dv": vd_items,
            "run_multinomial": bool(force_multinomial),
            "observation_cycle": observation_cycle,
            "output_dir": str(self.last_output_dir) if self.last_output_dir else "",
            "requested_at": datetime.now().isoformat(timespec="seconds"),
            "run_timestamp": run_timestamp,
            "program_name": PROGRAM_NAME,
            "computer_id": socket.gethostname(),
            "custom_formula": self.formula_text.get("1.0", "end").strip() or None,
        }

    def _run_r_subprocess(self, config_path: Path, output_dir: Path) -> None:
        rscript = self._resolve_rscript_path()
        run_multinomial = False
        try:
            cfg = json.loads(config_path.read_text(encoding="utf-8"))
            run_multinomial = bool(cfg.get("run_multinomial", False))
        except Exception:
            run_multinomial = False

        cmd = [
            rscript,
            str(RUN_ANALYSIS_SCRIPT),
            "--config",
            str(config_path),
            "--run_multinomial",
            "TRUE" if run_multinomial else "FALSE",
        ]

        self._log(f"Command: {' '.join(cmd)}")

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                cwd=str(APP_DIR),
            )

            if proc.stdout.strip():
                self._log("[R stdout]\n" + proc.stdout.strip())
            if proc.stderr.strip():
                self._log("[R stderr]\n" + proc.stderr.strip())

            if proc.returncode != 0:
                self._set_busy(False)
                self._show_error_async(
                    "R analysis failed",
                    f"Rscript exited with code {proc.returncode}.\n"
                    f"Check the execution log for details.",
                )
                return

            summary_file = output_dir / "results_summary.json"
            if summary_file.exists():
                self._log(f"Summary: {summary_file}")

            self._set_busy(False)
            self._show_info_async("Done", f"Analysis complete.\nOutput folder:\n{output_dir}")

        except subprocess.TimeoutExpired:
            self._set_busy(False)
            self._show_error_async(
                "Timeout",
                "R analysis exceeded 5 minutes. Check data/model complexity and try again.",
            )
        except FileNotFoundError:
            self._set_busy(False)
            self._show_error_async(
                "Rscript not found",
                "Unable to find Rscript.exe. Use 'Set Rscript Path' and retry.",
            )
        except Exception as exc:
            self._set_busy(False)
            self._show_error_async("Unexpected error", str(exc))

    def _open_last_output(self) -> None:
        if not self.last_output_dir or not self.last_output_dir.exists():
            messagebox.showinfo("No output", "No output directory is available yet.")
            return

        report_pdf = next(self.last_output_dir.glob("report_*.pdf"), None)
        report_html = next(self.last_output_dir.glob("report_*.html"), None)

        target = report_pdf or report_html or self.last_output_dir

        try:
            subprocess.Popen(["explorer", str(target)])
        except Exception as exc:
            messagebox.showerror("Error", f"Unable to open output:\n{exc}")

    def _resolve_rscript_path(self) -> str:
        cfg = self._load_config()
        explicit = cfg.get("rscript_path", "").strip()
        if explicit:
            return explicit
        return "Rscript"

    def _load_config(self) -> dict:
        if not CONFIG_PATH.exists():
            return {}
        try:
            return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _save_config(self, cfg: dict) -> None:
        CONFIG_PATH.write_text(json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")

    def _set_busy(self, busy: bool) -> None:
        def _apply() -> None:
            self.analyze_button["state"] = "disabled" if busy else "normal"
            self.status_var.set("Running analysis..." if busy else "Idle")

        self.root.after(0, _apply)

    def _show_error_async(self, title: str, msg: str) -> None:
        self.root.after(0, lambda: messagebox.showerror(title, msg))

    def _show_info_async(self, title: str, msg: str) -> None:
        self.root.after(0, lambda: messagebox.showinfo(title, msg))

    def _log(self, msg: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")

        def _append() -> None:
            self.log_text.insert("end", f"[{timestamp}] {msg}\n")
            self.log_text.see("end")

        self.root.after(0, _append)


def main() -> None:
    root = tk.Tk()
    app = RGLMMAnalyzerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()