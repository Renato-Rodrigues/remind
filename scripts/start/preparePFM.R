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

preparePFM <- function(cfg, verbose = TRUE) {
  say <- function(...) if (isTRUE(verbose)) cat("[preparePFM] ", ..., "
", sep = "")

  # Settings come from cfg$pfm (default.cfg, overridable per scenario). Every field
  # has a working default, so a coupled run needs NO environment set up at all. The
  # PFM_* variables still win if present - handy for a one-off experiment - but
  # nothing depends on them, and a fresh shell is never a reason for a run to fail.
  p <- cfg$pfm %||% list()
  gv <- function(field, env, default) {
    v <- Sys.getenv(env, "")
    if (nzchar(v)) return(v)
    if (!is.null(p[[field]]) && nzchar(as.character(p[[field]]))) return(p[[field]])
    default
  }
  sourceDir <- gv("source", "PFM_SOURCE", "../pfm-data")

  # --- settings DERIVED from the REMIND run, not restated -----------------------
  # Anything REMIND already knows is taken from REMIND. A second place to write the
  # same fact is a second place for it to be wrong, and this coupling has already
  # produced one silently-wrong configuration that way.
  #
  # SSP: cm_GDPpopScen is the run's GDP/population scenario, and the final-energy
  # weights must be projected under the SAME one or they describe a different world.
  ssp <- as.character(cfg$gms$cm_GDPpopScen %||% "")
  ssp <- if (nzchar(ssp)) sub("^gdp_", "", ssp) else gv("ssp", "PFM_SSP", "SSP2")
  #
  # Region mapping: the run's own, so the coupling can never deliver at a different
  # resolution than the model solves at.
  rmap <- basename(as.character(cfg$regionmapping %||% ""))
  if (!nzchar(rmap)) rmap <- gv("regionmapping", "PFM_MAPPING", "regionmapping_21_EU11.csv")
  #
  # Reference gdx: REMIND already copies path_gdx_ref into the run folder as
  # input_ref.gdx, so bind mode 2 needs no new switch and no copy of ours - it is
  # already there, already relative, already self-contained.
  refRun <- file.path(cfg$results_folder, "input_ref.gdx")

  coupled <- identical(as.character(cfg$gms$cm_taxCO2_regiDiff %||% ""), "11")
  if (!coupled) return(invisible(FALSE))

  dest <- file.path(cfg$results_folder, "pfm")

  # Which Run-Group? Normally there is exactly one prepared in pfm-data, so asking
  # the user to name it again is asking them to repeat a decision already made when
  # the folder was assembled - and to keep it in sync forever after. Auto-detect, and
  # only demand an answer when the folder is genuinely ambiguous.
  say("looking for PFM data under '", normalizePath(sourceDir, mustWork = FALSE),
      "' (cfg$pfm$source = '", sourceDir, "', wd = '", getwd(), "')")

  # Where the Run-Group lives. Two layouts, in order of preference:
  #   1. <source>/            - the files sitting directly in pfm-data (simplest)
  #   2. <source>/<group>/    - one or more Run-Group directories side by side
  # A Run-Group is identified by selected-models-psm.yml, not by its name, so neither
  # layout needs the group spelled out anywhere unless there are several to choose from.
  marker <- "selected-models-psm.yml"
  group <- gv("group", "PFM_GROUP", "")

  if (nzchar(group)) {
    src <- file.path(sourceDir, group)
    if (!file.exists(file.path(src, marker))) {
      stop("preparePFM: cfg$pfm$group = '", group, "' but no ", marker,
           " at '", src, "'.")
    }
  } else if (file.exists(file.path(sourceDir, marker))) {
    src <- sourceDir
    group <- basename(normalizePath(sourceDir, mustWork = FALSE))
    say("layout: files directly in the source folder")
  } else {
    cand <- list.dirs(sourceDir, full.names = FALSE, recursive = FALSE)
    cand <- cand[vapply(cand, function(g)
      file.exists(file.path(sourceDir, g, marker)), logical(1))]
    if (length(cand) == 1L) {
      group <- cand
      src <- file.path(sourceDir, group)
      say("layout: Run-Group auto-detected: ", group)
    } else if (length(cand) > 1L) {
      stop("preparePFM: ", length(cand), " Run-Groups in '", sourceDir, "' (",
           paste(cand, collapse = ", "), "). Set cfg$pfm$group to choose.")
    } else {
      stop("preparePFM: no PFM Run-Group found. Looked for '", marker, "' in '",
           sourceDir, "' and in each directory directly under it. ",
           "cfg$pfm$source is currently '", sourceDir, "'.")
    }
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
  # Only now, once every source file is confirmed present. Creating the folder earlier
  # left an EMPTY pfm/ behind on failure, which reads as "the copy worked" - the most
  # misleading possible state to leave a run folder in.
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(src, need), file.path(dest, need), overwrite = TRUE)

  # The panel the deployed spec was fitted on, addressed by the hash in the manifest.
  mf <- jsonlite::read_json(file.path(dest, "manifest.json"))
  panel <- paste0("panel_", mf$panel_hash, ".rds")
  # panels/ sits beside the group in the nested layout, and inside it in the flat one.
  # Exactly ONE panel is ever needed - the one named by panel_hash in manifest.json -
  # so a panels/ folder holding a single file is optional ceremony. Accept it flat
  # beside the other artifacts first, then the panels/ subfolder for a store copied
  # straight from a pfm Run-Group.
  panelCand <- c(file.path(src, panel),
                 file.path(sourceDir, panel),
                 file.path(src, "panels", panel),
                 file.path(sourceDir, "panels", panel))
  panelSrc <- panelCand[file.exists(panelCand)][1]
  if (is.na(panelSrc)) {
    stop("preparePFM: panel '", panel, "' not found in any of: ",
         paste(panelCand, collapse = ", "),
         " - the Run-Group and the panel store are out of sync.")
  }
  # In the RUN folder it always lands in pfm/panels/, because that is where
  # iterativePFM(modelDir = "pfm") looks. Only the SOURCE layout is flexible.
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
    paste0("couplingMapping: ", rmap),
    paste0("gdxRegionMapping: ", rmap),
    paste0("weightScenario: ", ssp),
    # weightYear: the year whose country-size distribution sets the within-region
    # aggregation weights. Kept NEAR-TERM on purpose - see default.cfg.
    paste0("weightYear: ", gv("weightYear", "PFM_WEIGHT_YEAR", 2025))
  )
  if (file.exists(refRun)) {
    cfgLines <- c(cfgLines, "refGdx: input_ref.gdx")
    say("reference price path: input_ref.gdx (from path_gdx_ref)")
  } else {
    say("NOTE: no input_ref.gdx - bind mode 2 will refuse to run without it")
  }
  say("derived from the run: ssp = ", ssp, ", regionmapping = ", rmap)
  writeLines(cfgLines, file.path(cfg$results_folder, "pfm-coupling.yml"))

  say("run folder is self-contained: ", dest)
  say("copied ", length(need), " artifacts + ", panel)
  invisible(TRUE)
}

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
