# FoodWebs.jl

Structural food web models, random community matrices, stability metrics and elliptic-law
predictions. The research code in
[ecological-reactivity](https://github.com/gabrielgellner/ecological-reactivity) builds on it.

## Layout

- `src/topology.jl`: `PredationMatrix` and the network generators (food chain, May, cascade,
  generalized cascade, niche, generalized niche), with Williams & Martinez style acceptance.
- `src/utils.jl`: connectance, basal and top species, trophic levels,
  `disconnected_from_basal`.
- `src/properties.jl`: the structural properties of Williams & Martinez (2000) and the
  characteristic path length, tested against their published tables.
- `src/biomass.jl`: observed biomass from node traits, steady-state biomass from energy flow,
  and pyramid measures.
- `src/feasible.jl`: `FeasibleCommunities`, the polytope of Lotka–Volterra parameters with a
  given equilibrium, its structure without a solver, and priors on it.
- `src/compartments.jl`: modularity of partitions of webs and weight matrices, and partition
  search (Louvain, simulated annealing, consensus), modularity z-scores against rewired and
  structural null webs, species roles and compartment homogeneity; `test/fixtures/make_modularity.py` checks
  modularity against NetworkX and igraph, and `make_modularity_optimum.py` gives exact optima.
- `src/motifs.jl`: three- and four-species motif census (mfinder IDs), species roles, and the
  rewiring null model for z-scores; `test/fixtures/make_motif_census.py` checks the census
  against igraph and NetworkX.
- `src/data.jl`: the `FoodWeb` type and readers for the empirical webs in `data/` (see
  `data/README.md`).
- `src/randcmat.jl`: random community matrices with a web's structure (`random_predprey`).
- `src/stability.jl`: spectral and numerical abscissa.
- `src/elliptic_law.jl`: elliptic-law predictions and off-diagonal pair utilities.
- `test/`: one test file per source file, each testset with its own seeded `Xoshiro`.
- `docs/`: the Documenter site (own environment), deployed to GitHub Pages by
  `.github/workflows/Documentation.yml` on pushes to main and release tags.

## Conventions

- Food webs: `web[i, j] == 1` means species `j` eats species `i` (prey are rows, predators
  are columns). `PredationMatrix` stores a `Matrix{Bool}`.
- Connectance is `L / S^2`.
- Community matrices use the Jacobian convention: `J[i, j]` is the effect of species `j` on
  the growth of species `i`, so for predator `j` and prey `i`, `J[i, j] < 0 < J[j, i]`.
- Random generators take an optional leading `rng::AbstractRNG`, and accept a web only when
  it is connected and its connectance is within `tol` of the target (`basal_paths = true`
  also requires every species to have a feeding path up from a basal species).
- Julia 1.13 only. Imports are explicit (`using X: a, b`), checked by ExplicitImports.

## Tests

Test model properties rather than single outputs: connectance unbiased with acceptance off,
structural invariants (acyclic cascades, interval diets), sign structure, reproducibility
from a seeded rng, argument errors, hand-derived values. Check a new stochastic test under
several shifted seed sets before trusting it; several early tests here passed only for the
seed they were written with.

## Commands

```sh
just setup      # install Runic and JETLS as Julia apps, instantiate environments, enable hooks
just fmt        # format with Runic
just ci         # formatting check, JETLS lint, tests (run by the pre-push hook)
just gallery    # render the Quarto notes in gallery/ (webs, spectra, equilibria)
just docs       # build the Documenter site to docs/build
```

`main` is protected: changes go through a PR. GitHub Actions runs nothing automatically; run
`just ci` before pushing or merging. The pre-push hook (`.githooks/pre-push`, enabled by
`just setup`) runs it on every branch push. The CI workflow can still be run by hand.

## Commits

Conventional commits, `type(scope): subject`, with the issue as `(#N)` at the end when
there is one. git-cliff builds `CHANGELOG.md` from the subjects (`cliff.toml`): the type
picks the section an entry lands in, and git-cliff drops a subject it cannot parse.

- `feat` → Features, `fix` → Bug Fixes, `docs` → Documentation, `perf` → Performance,
  `refactor` → Refactoring, `test` → Testing, `chore` → Miscellaneous.
- `!` after the type or scope (`fix!:`, `feat(topology)!:`) marks a breaking change: it heads
  the entry under Breaking Changes and bumps the minor version while we are pre-1.0.
- `ci`, `chore(release)` and merge commits are left out of the changelog.
- Write the subject as the changelog line a user reads: what changed, in the imperative
  ("Add the nested-hierarchy model").

Releases run through `just release-1-preview` … `just release-4-merge`: preview the
version, write the changelog on a `release/vX` branch (edit it by hand if needed), open the
release PR with the version bumped in `Project.toml`, then merge, tag the merge commit and
create the GitHub release from the edited changelog section.

## Writing style

These rules cover all prose written for this repo: comments, docstrings, documentation,
commit messages, issues and PRs.

Write in plain language. State claims directly, in the order the reader
needs them. Use subject-verb-object sentences, concrete nouns, and the
actual number, name, or mechanism instead of an allusion to it.

Follow these rules from Strunk & White:

- Put statements in positive form.
- Use definite, specific, concrete language.
- Omit needless words.

Avoid these constructions everywhere (sentences, headings, titles, lists):

- Contrast by negation: "not X, but Y", "it's not about X, it's about Y",
  "X isn't the problem; Y is". State Y. Mention X only if someone
  actually holds that view.
- Hook clauses that read like half a question: "what X doesn't tell you",
  "the part nobody mentions", "the evidence Y lacks", "why this matters".
  State the fact the clause points at.
- Colon reveals: "The real issue: ...", "Here's the thing: ...",
  "The short version: ...". Drop the setup and write the sentence.
- Punchline structure: ending a sentence or list on a twist, a dramatic
  final item, or a short fragment for emphasis.
- Rhetorical questions that the text then answers.
- Framing work or ideas by what is missing when you can name what is
  present or proposed.

Titles (milestones, issues, PRs, headings):

- Name the deliverable or scope as a plain noun phrase or imperative.
  A reader should know from the title alone what gets merged or decided.
- Use a colon only to separate a component name from the work item
  (e.g. "Query DSL: add DuckDB backend"), never to introduce a hook.

Examples:

```text
Bad:  Conformance evidence: what the green checks do not prove
Good: Conformance: add negative and edge-case test fixtures

Bad:  Query DSL: a second backend, and the evidence the first one lacks
Good: Query DSL: add second backend; add validation tests for SQL backend

Bad:  The green checks tell you less than you'd think.
Good: The green checks only cover valid inputs; there are no tests
      for malformed DSDs or constraint violations.

Bad:  The problem isn't performance — it's correctness.
Good: The cache returns stale codelists after a version bump.

Bad:  Here's what the benchmark doesn't show: cold-start latency.
Good: The benchmark measures warm queries only; cold-start latency
      is untested.
```

Before finishing, reread the output. Rewrite any sentence that delays its
point, contrasts with a claim nobody made, or alludes to a fact instead of
stating it.
