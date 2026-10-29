# Repository workflow

- Work in the `franceslinyc/mtd` fork. Do not modify, push to, or open PRs against `xzheng42/mtd` unless Frances explicitly asks.
- Before editing, inspect the current branch, working tree, target files, and relevant PR base. Preserve unrelated user edits and generated files.
- The main project goal is to merge the user's source implementation into this fork. Stay as close to that code as possible; omit comments, but preserve required third-party attribution. Do not change implementation code to "fix" suspected errors without describing the issue and getting Frances's approval first.
- Add or update tests with behavior changes. Run the relevant tests before reporting completion; use `SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.sdk` for R package builds on this Mac.
- For each focused change, create one new branch from current `main` and one PR targeting `main`; avoid stacked or multiple overlapping PRs. Wait for Frances to review, approve, and merge the PR. Only after it is merged, sync local `main` and delete that PR's local and remote feature branches.
- Do not merge PRs before Frances approves, delete branches before their PR is merged, discard user changes, or clean generated files unless Frances explicitly asks.
- Before creating a PR, summarize its scope, test results, and any known limitations. Leave merge decisions to Frances unless explicitly instructed otherwise.
