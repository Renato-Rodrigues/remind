# Make a REMIND run folder self-contained for the PFM coupling.
#
# The coupling is invoked from GAMS as `Rscript -e 'pfm::iterativePFM()'`, whose
# working directory IS the run folder. Everything it needs therefore lives under
# <run>/pfm/ and is addressed RELATIVELY, so the folder can be moved, archived or
# re-run on another machine without editing a path. Nothing outside the run folder is
# read at coupling time.
#
# This replaces putting PFM settings in REMIND's own .Rprofile: that file is shared
# with every non-coupled run, is not per-scenario, and would have carried absolute
# paths. Called from submit.R after the standard files2export copy.

preparePFM <- function(cfg, sourceDir = Sys.getenv("PFM_SOURCE", "../pfm-data"),
                       verbose = TRUE) {
  say <- function(...) if (isTRUE(verbose)) cat("[preparePFM] ", ..., "\n", sep = "")

  coupled <- identical(as.character(cfg$gms$cm_taxCO2_regiDiff %||% ""), "11")
  if (!coupled) return(invisible(FALSE))

  dest <- file.path(cfg$results_folder, "pfm")
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)

  group <- Sys.getenv("PFM_GROUP", "psm-country-v3")
  src <- file.path(sourceDir, "output", group)
  if (!dir.exists(src)) {
    stop("preparePFM: PFM Run-Group not found at '", src, "'. Set PFM_SOURCE to the ",
         "directory holding output/<group>/, or copy the group in by hand. The run ",
         "cannot be made self-contained without it.")
  }

  # Only what iterativePFM() actually reads. Copying the whole Run-Group would drag in
  # sweep.rds and the projection fan-out - hundreds of MB of things the coupling never
  # opens - into every run folder.
  need <- c("selected-models-psm.yml", "manifest.json", "frontier.rds",
            "temporal-validation.rds",
            "donor-assignment-band-Bulk.rds", "donor-assignment-band-Diffuse.rds")
  missing <- need[!file.exists(file.path(src, need))]
  if (length(missing)) {
    stop("preparePFM: the Run-Group is missing ", paste(missing, collapse = ", "),
         ". The band assignments come from analysis/psm-donor-assumptions.R; without ",
         "them the coupling refuses to run rather than reverting to phi = 1.")
  }
  file.copy(file.path(src, need), file.path(dest, need), overwrite = TRUE)

  # The panel the deployed spec was fitted on, addressed by the hash in the manifest.
  mf <- jsonlite::read_json(file.path(dest, "manifest.json"))
  panel <- paste0("panel_", mf$panel_hash, ".rds")
  panelSrc <- file.path(sourceDir, "output", "panels", panel)
  if (!file.exists(panelSrc)) {
    stop("preparePFM: panel '", panel, "' not found at ", panelSrc,
         " - the Run-Group and the panel store are out of sync.")
  }
  dir.create(file.path(dest, "panels"), showWarnings = FALSE)
  file.copy(panelSrc, file.path(dest, "panels", panel), overwrite = TRUE)

  # Static settings, RELATIVE to the run folder. bindMode and theta are deliberately
  # absent: presolve.gms writes them to pfm-coupling-runtime.yml from the scenario
  # config, so there is exactly one place they are set.
  cfgLines <- c(
    "# Written by scripts/start/preparePFM.R - do not edit by hand.",
    "# All paths are RELATIVE to the run folder, which is the working directory when",
    "# GAMS calls Rscript. bindMode and theta come from pfm-coupling-runtime.yml.",
    paste0("group: ", group),
    "resultsDir: pfm",
    "modelDir: pfm",
    "couplingMapping: regionmapping_21_EU11.csv",
    "gdxRegionMapping: regionmapping_21_EU11.csv",
    paste0("weightScenario: ", Sys.getenv("PFM_SSP", "SSP2")),
    "weightYear: 2050"
  )
  refGdx <- Sys.getenv("PFM_REF_GDX", "")
  if (nzchar(refGdx)) {
    # Copy the reference gdx IN, so bind mode 2 does not depend on an outside path.
    file.copy(refGdx, file.path(dest, "reference.gdx"), overwrite = TRUE)
    cfgLines <- c(cfgLines, "refGdx: pfm/reference.gdx")
    say("reference gdx copied in for bind mode 2")
  }
  writeLines(cfgLines, file.path(cfg$results_folder, "pfm-coupling.yml"))

  say("run folder is self-contained: ", dest)
  say("copied ", length(need), " artifacts + ", panel)
  invisible(TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
