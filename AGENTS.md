# Repository workflow

- Work in the `franceslinyc/mtd` fork. Do not modify, push to, or open PRs against `xzheng42/mtd` unless Frances explicitly asks.
- Before editing, inspect the current branch, working tree, target files, and relevant PR base. Preserve unrelated user edits and generated files.
- Keep each change focused. When merging code from `code_mtd_v17_b4`, review one file or coherent chunk at a time; omit source comments unless specifically requested, while retaining necessary attribution.
- Add or update tests with behavior changes. Run the relevant tests before reporting completion; use `SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk` for R package builds on this Mac.
- Prefer separate PRs based on `main`; avoid stacked PRs unless a dependency makes stacking necessary and Frances agrees. If stacking, state each PR's base branch clearly and verify where commits land; a PR merged into a feature branch is not yet on `main`.
- Do not merge PRs, delete local or remote branches, discard user changes, or clean generated files unless Frances explicitly asks.
- Before creating a PR, summarize its scope, test results, and any known limitations. Leave merge decisions to Frances unless explicitly instructed otherwise.
