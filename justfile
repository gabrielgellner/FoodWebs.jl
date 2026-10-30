# Development tasks for FoodWebs.jl. Run `just` to list recipes.

# Julia apps (runic, jetls) are installed to ~/.julia/bin
export PATH := home_directory() / ".julia" / "bin" + ":" + env_var("PATH")

julia := "julia --startup-file=no"

# List recipes
default:
    @just --list

# Install the formatter and linter as Julia apps and instantiate all environments
setup:
    {{ julia }} -e 'using Pkg; Pkg.Apps.add("Runic"); Pkg.Apps.add(url = "https://github.com/aviatesk/JETLS.jl", rev = "release")'
    {{ julia }} --project -e 'using Pkg; Pkg.instantiate()'

# Format all Julia files in place with Runic
fmt:
    runic --inplace .

# Check formatting without changing files, printing a diff
fmt-check:
    runic --check --diff .

# Static analysis with JETLS (JET inference plus lowering lints)
lint:
    jetls check -q .

# Run the test suite (Aqua, ExplicitImports and unit tests)
test:
    {{ julia }} --project -e 'using Pkg; Pkg.test()'

# Run every check that CI runs
check: fmt-check lint test
