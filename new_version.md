# Releasing a New Version

## 1. Open the PR

```sh
gh pr create \
  --base main \
  --title "vX.Y.Z: short description" \
  --body "$(cat <<'EOF'
## Summary
- Change 1
- Change 2
EOF
)"
```

## 2. Watch CI and merge

Open the PR in the browser:

```sh
gh pr view --web
```

Wait for all CI jobs to pass, then squash-merge and delete the branch.

## 3. Update local main

```sh
git checkout main
git pull
```

## 4. Register in the private registry

From the Julia REPL:

```julia
using LocalRegistry
register("TransitionFunctionIteration", registry="MacroFiscalAPRegistry", push=true)
```

This adds the new version to `~/.julia/registries/MacroFiscalAPRegistry/` and pushes to `git@github.com:elenev/MacroFiscalAP-registry.git`. Requires SSH key access.

## 5. TagBot tags automatically

TagBot runs daily (`0 0 * * *`) and creates a tag and GitHub release once it detects the new version in the registry. To trigger immediately:

```sh
gh workflow run TagBot.yml
```
