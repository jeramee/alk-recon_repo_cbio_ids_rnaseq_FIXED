import json
import subprocess
from collections import defaultdict
from pathlib import Path

# VS Code can dump diagnostics via "code --status" doesn't include them,
# so we read the Problems view by asking VS Code to export via "Developer: Open Diagnostics"
# ...but VS Code doesn't expose a direct CLI for that.
# The reliable way is: use the Problems panel "Copy All" (Option 1).

print("VS Code does not expose Problems diagnostics via CLI in a stable way.")
print("Use Option 1 (Problems panel -> ... -> Copy All) for the raw dump.")
