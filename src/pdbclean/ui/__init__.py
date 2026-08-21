"""Web UI for the COMP702 PDBClean pipeline.

The UI is a thin presentation layer over the same backend the CLI uses:
:mod:`pdbclean.runconfig` resolves configuration, :mod:`pdbclean.pipeline`
plans and gates stages, :mod:`pdbclean.duplicates` queries results and
:mod:`pdbclean.run_provenance` records runs.

No scientific logic lives in this package or in its JavaScript.
"""

from pdbclean.ui.server import serve

__all__ = ["serve"]
