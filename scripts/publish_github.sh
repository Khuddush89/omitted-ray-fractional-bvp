#!/bin/bash
# Creates a PUBLIC repository; never force-pushes or overwrites remote history.
set -euo pipefail
cd "$(dirname "$0")/.."
repo="Khuddush89/omitted-ray-fractional-bvp"
url="https://github.com/$repo.git"
command -v git >/dev/null || { echo "Install Apple command-line tools: xcode-select --install"; exit 1; }
command -v gh >/dev/null || { echo "Install GitHub CLI: brew install gh (or visit https://cli.github.com/)"; exit 1; }
if ! gh auth status --hostname github.com >/dev/null 2>&1; then
  gh auth login --hostname github.com --git-protocol https --web
fi
login="$(gh api --hostname github.com user --jq .login)"
if [ "$login" != "Khuddush89" ]; then
  echo "Please sign in as Khuddush89. Active account: $login"
  echo "Run: gh auth switch --hostname github.com --user Khuddush89"
  exit 1
fi
gh auth setup-git --hostname github.com
if [ ! -d .git ]; then
  if gh repo view "$repo" >/dev/null 2>&1; then
    echo "The remote repository already exists. Stopping to preserve its contents."
    echo "Use its existing checkout to review and apply this package."
    exit 1
  fi
  git init
  git symbolic-ref HEAD refs/heads/main
fi
branch="$(git symbolic-ref --short HEAD)"
[ "$branch" = main ] || { echo "Expected branch main; found $branch. Stopping."; exit 1; }
if git remote get-url origin >/dev/null 2>&1; then
  origin="$(git remote get-url origin)"
  case "$origin" in
    "$url"|"https://github.com/$repo"|"git@github.com:$repo.git") ;;
    *) echo "Unexpected origin: $origin. Stopping."; exit 1 ;;
  esac
fi
if ! git config user.name >/dev/null; then git config user.name 'Mahammad Khuddush'; fi
if ! git config user.email >/dev/null; then
  uid="$(gh api --hostname github.com user --jq .id)"
  git config user.email "$uid+$login@users.noreply.github.com"
fi
git add README.md LICENSE CITATION.cff requirements.txt .gitignore src scripts results figures assets UPLOAD_TO_GITHUB.md
if ! git diff --cached --quiet; then
  git commit -m "Present fractional BVP code, data and illustrated research guide"
fi
if ! git remote get-url origin >/dev/null 2>&1; then
  gh repo create "$repo" --public --source=. --remote=origin \
    --description "Fractional BVPs: Green functions, positive solutions, omitted ray fixed point theorem, Python code and figures."
fi
git push -u origin main
if ! gh repo edit "$repo" --add-topic fractional-calculus --add-topic boundary-value-problems \
  --add-topic green-functions --add-topic fixed-point-theory --add-topic numerical-analysis --add-topic python; then
  echo "Upload succeeded; topics could not be added automatically."
fi
echo "Published: https://github.com/$repo"
if command -v open >/dev/null; then open "https://github.com/$repo"; fi
