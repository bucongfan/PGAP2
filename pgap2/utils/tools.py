"""
Backward-compatible shim for pgap2.utils.tools.

This module has been renamed to pgap2.utils.graph_utils.
All imports are re-exported here for backward compatibility.
New code should import from pgap2.utils.graph_utils directly.
"""
from pgap2.utils.graph_utils import *  # noqa: F401, F403
