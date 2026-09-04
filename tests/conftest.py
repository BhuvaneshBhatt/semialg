from __future__ import annotations

import os
import random


def pytest_collection_modifyitems(session, config, items):
    """Optionally randomize test execution order with a reproducible seed.

    Normal local/default runs retain source order. CI can set
    ``SEMIALG_RANDOM_TEST_ORDER`` to expose process-state and cache-order bugs.
    """

    seed = os.environ.get("SEMIALG_RANDOM_TEST_ORDER")
    if seed is None:
        return
    rng = random.Random(int(seed))
    rng.shuffle(items)
