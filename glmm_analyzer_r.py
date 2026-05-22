"""
R-backed GLMM analyzer GUI (first robust version).

This app keeps the familiar Tkinter workflow but delegates model fitting,
diagnostics, emmeans, and report generation to R via Rscript.


24 avril 2026, add save as the input file name instead of timestamp to make it easier to track outputs from different datasets.
"""

from __future__ import annotations

import json
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

IV_ROLE_VALUES = ["fixed_effect", "random_intercept", "random_effect", "random_slope"]
IV_SCOPE_VALUES = ["within_subject", "between_subject"]
DV_TYPE_VALUES = ["counts", "rt", "proportion", "proportion_without_non_visible", "proportion_multinomial"]


class RGLMMAnalyzerGUI:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("GLMM Analyzer (R Backend)")
        self.root.geometry("1200x920")

        self.data_path: Path | None = None
        self.data_df: pd.DataFrame | None = None
        self.sheet_name: str | None = None
        self.last_output_dir: Path | None = None

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

        for i in range(4):
            row = ttk.Frame(vi_box)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"VI {i + 1}:", width=8).pack(side="left")
            vi_combo = ttk.Combobox(row, state="readonly", width=28)
            vi_combo.pack(side="left", padx=4)
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
            self.vi_role_combos.append(vi_role)

            ttk.Label(row, text="Scope:", width=7).pack(side="left", padx=(10, 2))
            vi_scope = ttk.Combobox(
                row,
                state="readonly",
                width=16,
                values=IV_SCOPE_VALUES,
            )
            vi_scope.set("within_subject")
            vi_scope.pack(side="left")
            self.vi_scope_combos.append(vi_scope)

        # VDs + type (+ non-visible for dedicated proportion mode)
        vd_box = ttk.LabelFrame(frame, text="Dependent Variables (up to 10)", padding=6)
        vd_box.pack(fill="x", pady=6)
        self.vd_combos: list[ttk.Combobox] = []
        self.vd_type_combos: list[ttk.Combobox] = []

        cycle_row = ttk.Frame(vd_box)
        cycle_row.pack(fill="x", pady=(0, 4))
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
        non_visible_row.pack(fill="x", pady=(0, 4))
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
            self.vd_combos.append(vd_combo)

            ttk.Label(row, text="Type:", width=6).pack(side="left", padx=(10, 2))
            vd_type = ttk.Combobox(
                row,
                state="readonly",
                width=30,
                values=DV_TYPE_VALUES,
            )
            vd_type.set("counts")
            vd_type.pack(side="left")
            self.vd_type_combos.append(vd_type)

            vd_type.bind("<<ComboboxSelected>>", lambda _e: self._on_vd_type_changed())

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

    def _load_columns_for_selected_sheet(self) -> None:
        if not self.data_path:
            return

        sheet = self.sheet_combo.get().strip()
        if not sheet:
            return

        try:
            self.sheet_name = sheet
            self.data_df = pd.read_excel(self.data_path, sheet_name=sheet, engine="openpyxl")
            columns = [""] + [str(c) for c in self.data_df.columns]

            self.id_combo["values"] = columns
            for combo in self.vi_combos:
                combo["values"] = columns
            for combo in self.vd_combos:
                combo["values"] = columns
            self.global_nonvisible_combo["values"] = columns
            self._update_global_nonvisible_state()

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
            combo.get().strip().lower() in ("proportion_without_non_visible", "proportion_multinomial")
            for combo in self.vd_type_combos
        )
        if requires_non_visible:
            self.global_nonvisible_combo.configure(state="readonly")
        else:
            self.global_nonvisible_combo.set("")
            self.global_nonvisible_combo.configure(state="disabled")

    def _on_vd_type_changed(self) -> None:
        self._update_global_nonvisible_state()

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

    def _collect_payload(self) -> dict:
        id_col = self.id_combo.get().strip()

        try:
            observation_cycle = int(self.obs_cycle_var.get().strip())
        except ValueError as exc:
            raise ValueError("Observation cycle must be an integer (e.g., 30).") from exc

        if observation_cycle <= 0:
            raise ValueError("Observation cycle must be > 0.")

        vi_items: list[dict] = []
        for vi_combo, vi_role_combo, vi_scope_combo in zip(
            self.vi_combos,
            self.vi_role_combos,
            self.vi_scope_combos,
        ):
            vi = vi_combo.get().strip()
            if vi:
                scope = vi_scope_combo.get().strip().lower() or "within_subject"
                vi_items.append(
                    {
                        "name": vi,
                        "model_role": vi_role_combo.get().strip().lower() or "fixed_effect",
                        "subject_scope": scope,
                        # Backward-compatible field consumed by older R code paths.
                        "type": "within" if scope == "within_subject" else "between",
                    }
                )

        if not vi_items:
            raise ValueError("Please select at least one independent variable.")

        vd_items: list[dict] = []
        force_multinomial = False
        shared_non_visible_col = self.global_nonvisible_combo.get().strip()
        for vd_combo, vd_type_combo in zip(
            self.vd_combos,
            self.vd_type_combos,
        ):
            vd = vd_combo.get().strip()
            if vd:
                dv_type = vd_type_combo.get().strip().lower()
                if dv_type == "proportion_multinomial":
                    # Backward-compatible: reuse existing proportion-without-non-visible pipeline.
                    # Multinomial section is toggled globally via run_multinomial.
                    force_multinomial = True
                    dv_type = "proportion_without_non_visible"
                item: dict[str, str] = {"name": vd, "type": dv_type}
                if dv_type == "proportion_without_non_visible":
                    if not shared_non_visible_col:
                        raise ValueError(
                            "Please select the shared Non-visible column for DVs using proportion_without_non_visible/proportion_multinomial."
                        )
                    item["non_visible_column"] = shared_non_visible_col

                vd_items.append(item)

        if not vd_items:
            raise ValueError("Please select at least one dependent variable.")

        return {
            "input": str(self.data_path),
            "sheet": self.sheet_name,
            "animal_id": id_col,
            "iv": vi_items,
            "dv": vd_items,
            "run_multinomial": bool(force_multinomial),
            "observation_cycle": observation_cycle,
            "output_dir": str(self.last_output_dir) if self.last_output_dir else "",
            "requested_at": datetime.now().isoformat(timespec="seconds"),
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
