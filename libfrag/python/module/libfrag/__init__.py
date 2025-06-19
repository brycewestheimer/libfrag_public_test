"""
libfrag: A quantum chemistry library for fragment-based calculations
"""

try:
    from ._libfrag import *
except ImportError:
    # During development, the C++ extension might not be built yet
    print("Warning: C++ extension not found. Build the project first.")
    pass

__version__ = "0.1.0"

def pure_python():
    """
    Example pure Python function for testing
    """
    return "Hello from libfrag Python API!"