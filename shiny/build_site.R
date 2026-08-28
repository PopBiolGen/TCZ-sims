# build_site.R
#
# Export the Shiny app to a self-contained shinylive site under shiny/site/.
# The result is plain static files (HTML + JS + wasm) that run the app entirely
# in the visitor's browser via webR -- no R process, no Shiny Server. Serve it
# from any static host (Netlify / Cloudflare Pages / GitHub Pages / S3 / nginx).
#
# Run from the repo root, after prep_app_data.R has regenerated the app data:
#   Rscript shiny/prep_app_data.R      # writes shiny/data/app_data.rds
#   Rscript shiny/build_site.R         # writes shiny/site/
#
# The first run downloads ~400 MB of shinylive web assets (webR runtime + wasm
# package binaries) into a user-level cache (see shinylive::assets_info()).
# Subsequent runs reuse the cache and take a few seconds.

if (!requireNamespace("shinylive", quietly = TRUE)) {
  install.packages("shinylive", repos = "https://cloud.r-project.org")
}

app_r    <- "shiny/app.R"
app_data <- "shiny/data/app_data.rds"
out_dir  <- "shiny/site"

stopifnot(
  "shiny/app.R not found -- run from the repo root"       = file.exists(app_r),
  "shiny/data/app_data.rds not found -- run prep_app_data.R first" = file.exists(app_data)
)

# Export from a clean staging copy so shiny/ helper scripts (this file,
# prep_app_data.R, rsconnect/, the .Rproj) don't get bundled into the app.
staging <- tempfile("tcz-shinylive-")
dir.create(file.path(staging, "data"), recursive = TRUE)
file.copy(app_r, file.path(staging, "app.R"))
file.copy(app_data, file.path(staging, "data", "app_data.rds"))
on.exit(unlink(staging, recursive = TRUE), add = TRUE)

unlink(out_dir, recursive = TRUE)
shinylive::export(staging, out_dir)

cat(
  "\nExported to ", out_dir, "\n",
  "Preview locally:\n",
  "  Rscript -e 'httpuv::runStaticServer(\"", out_dir, "\", port = 8008)'\n",
  "  then open http://localhost:8008\n",
  sep = ""
)
