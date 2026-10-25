# Security policy

## Reporting a vulnerability

Please report security problems privately, through GitHub's
[private vulnerability reporting](https://github.com/obormot65536/tessera-studio/security/advisories/new)
(the **Security** tab, then **Report a vulnerability**). Don't open a public issue.

Examples:
- A project file or an imported model that makes Tessera run commands or write outside its
  working folder.
- A way to make a workflow in this repository leak a token.

Tessera is maintained by one person. Expect an answer within two weeks. Fixes land on `main`
and in the next release.

## Not a security problem

A wrong or misleading result, such as a strength ratio, stress or deflection, is a correctness
bug. It matters just as much here: please open an ordinary issue, with the project file if you
can share it.

## Supported versions

Only the latest release and `main` get fixes.

## Contributors

Changes under `.claude/` (agents, hooks, settings) and `.github/` run on the machines of anyone
who opens the repository in Claude Code, or in CI. A PR touching them gets an especially careful
review.
