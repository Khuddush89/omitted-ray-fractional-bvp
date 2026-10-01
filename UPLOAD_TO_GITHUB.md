# ☁️ Upload from your Mac

This package is prepared for **Khuddush89/omitted-ray-fractional-bvp**. It has not been published by ChatGPT. Running the helper below creates a **public** repository and uploads the project, including saved results and figures.

1. Save the revised `omitted-ray-fractional-bvp.zip` in Downloads. Extract it into a new folder so it does not mix with an older copy.
2. Install GitHub CLI if necessary. With Homebrew installed, run `brew install gh`. Otherwise use the macOS installation instructions at https://cli.github.com/.
3. Open Terminal in the extracted `omitted-ray-fractional-bvp` folder and run:

```bash
bash scripts/publish_github.sh
```

The helper opens GitHub sign-in when needed. Sign in as **Khuddush89**. If macOS requests Apple command-line tools, complete that installation and run the helper again.

The helper checks the account, initializes a local Git repository, commits the package, creates the public GitHub repository, pushes `main`, adds topic tags and opens the repository. It uses your existing Git author settings when present; otherwise it configures the author's name and GitHub no-reply email only in this repository.

If the repository already exists and this folder is not its checkout, the helper stops to preserve the existing content. It never force-pushes. If upload is interrupted after repository creation, rerun it from the same folder.

## Later changes

From the same project folder:

```bash
git status
git add README.md assets figures scripts results src
git commit -m "Update documentation and research outputs"
git push
```

## Rebuild the visual gallery

Original plots: `python scripts/make_figures.py`.
README previews: `bash scripts/export_previews.sh` (requires Poppler).
Banner and section headings: editable SVG files in `assets/`.

GitHub CLI references: https://cli.github.com/manual/gh_repo_create and https://cli.github.com/manual/gh_auth_login.
