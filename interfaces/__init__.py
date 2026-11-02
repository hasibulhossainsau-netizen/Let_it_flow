"""Plugin interfaces for swappable external tools.

Implementations are selected by name through the ``tools`` section of
config/config.yaml. To add a tool, subclass ``OffTargetSearcher`` or
``StructurePredictor`` and register it in the dictionaries below.
"""

from .base import HIT_COLUMNS, OffTargetSearcher, StructurePredictor
from .blastn import BlastnSearcher
from .cas_offinder import CasOFFinderSearcher
from .rnafold import RNAfoldPredictor

OFF_TARGET_SEARCHERS = {
    CasOFFinderSearcher.name: CasOFFinderSearcher,
    BlastnSearcher.name: BlastnSearcher,
}
STRUCTURE_PREDICTORS = {
    RNAfoldPredictor.name: RNAfoldPredictor,
}


def _lookup(registry: dict, name: str, kind: str):
    try:
        return registry[name]
    except KeyError:
        raise ValueError(
            f"Unknown {kind} '{name}'. Available: {', '.join(sorted(registry))}"
        ) from None


def get_structure_predictor(config: dict) -> StructurePredictor:
    """Build the structure predictor selected in ``config['tools']``."""
    tools = config["tools"]
    cls = _lookup(STRUCTURE_PREDICTORS, tools["structure_predictor"],
                  "structure predictor")
    return cls(timeout_s=tools["rnafold_timeout_s"])


def get_off_target_searchers(config: dict) -> list:
    """Build every off-target searcher listed in ``config['tools']``."""
    tools, guide, blast = config["tools"], config["guide"], config["blast"]
    searchers = []
    for name in tools["off_target_searchers"]:
        cls = _lookup(OFF_TARGET_SEARCHERS, name, "off-target searcher")
        if cls is CasOFFinderSearcher:
            searchers.append(cls(device=tools["cas_offinder_device"],
                                 pam=guide["cas_offinder_pam"],
                                 spacer_length=guide["spacer_length"]))
        elif cls is BlastnSearcher:
            searchers.append(cls(task=blast["task"],
                                 word_size=blast["word_size"],
                                 perc_identity=blast["perc_identity"],
                                 threads=config["threads"]))
        else:
            searchers.append(cls())
    return searchers
