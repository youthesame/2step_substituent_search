## Common Commands

```bash
# Environment setup
source .venv/bin/activate

# Generate CRYSTAL input files
uv run python src/*.py
```

## Documentation Standards

- **Session Startup**: At the beginning of each Claude session, create a TODO section using the YAML format from the template to plan tasks
- **Progress Tracking**: Update TODO status throughout the session as tasks move from `pending` → `in_progress` → `completed`
- **Session Documentation**: After completing implementations, create concise logs in `_docs/` using the template at `_docs/templates/yyyy-mm-dd_feature-name.md`. Use `date "+%Y-%m-%d %H:%M:%S"` for timestamps.
- **History Review**: Review existing logs in `_docs/` at project startup to understand implementation history and design patterns.
