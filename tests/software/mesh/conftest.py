"""Shared fixtures of the mesh tests: ``official`` (the official IEA-15 windIO v2 topology) is defined once in
``test_windio2_topology.py`` and exposed here so the other windIO v2 test modules request it without importing it."""

from tests.software.mesh.test_windio2_topology import official  # noqa: F401
