# Git Cheat-Sheet (SpaX)

A quick reference for everyday git. Your repo lives at
**https://github.com/geck-000/SpaX** (private).

## The everyday cycle

```bash
git status               # See what changed (safe, run anytime)
git add .                # Stage all changes for the next snapshot
git commit -m "message"  # Save a snapshot locally
git push                 # Upload your commits to GitHub
```

> Tip: `git add .` stages everything. To stage one file: `git add Spatium_Standalone.py`

## Looking around

```bash
git log --oneline        # Compact history of commits
git log --oneline -5     # Just the last 5
git diff                 # What you changed but haven't staged yet
git diff --staged        # What you've staged but not committed
git show HEAD            # Details of the most recent commit
```

## Undoing things

```bash
git restore FILE             # Discard unstaged changes to a file (careful!)
git restore --staged FILE    # Unstage a file (keeps your edits)
git revert <commit>          # Make a new commit that undoes an old one (safe)
git reset --soft HEAD~1      # Undo the last commit, keep the changes staged
```

## Branches (parallel lines of work)

```bash
git branch                   # List branches; * marks the current one
git switch -c my-feature     # Create and move to a new branch
git switch main              # Go back to the main branch
git merge my-feature         # Merge another branch into the current one
```

## Versions / releases (tags)

A **tag** marks a specific commit as a named version, e.g. `v1.0.0`.

```bash
git tag                              # List all tags
git tag -a v1.1.0 -m "what's new"    # Create an annotated version tag
git push origin v1.1.0               # Push one tag to GitHub
git push --tags                      # Push all tags
```

Version numbers follow **MAJOR.MINOR.PATCH** (semantic versioning):
- **PATCH** (v1.0.1) — bug fixes, small tweaks
- **MINOR** (v1.1.0) — new features, backwards-compatible
- **MAJOR** (v2.0.0) — big or breaking changes

## Syncing with GitHub

```bash
git push                 # Send your commits up
git pull                 # Bring down changes from GitHub (e.g. another machine)
```

## If you get stuck

- `git status` almost always tells you what to do next.
- Nothing is truly lost once committed — ask Claude and it can usually recover it.
