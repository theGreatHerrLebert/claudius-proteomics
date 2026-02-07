"""
Engine Registry

Maps engine names to their job classes. Adding a new engine requires:
1. Create runner/engines/newengine_job.py extending EngineJob
2. Add to ENGINES dict below
3. Create scripts/engine_parsers/newengine_parser.py extending BaseParser
4. Add config section to config/config.yaml
"""

from runner.engines.fragpipe_job import FragPipeJob
from runner.engines.diann_job import DiannJob
from runner.engines.sage_job import SageJob

ENGINES = {
    "fragpipe": FragPipeJob,
    "diann": DiannJob,
    "sage": SageJob,
}

__all__ = ["ENGINES", "FragPipeJob", "DiannJob", "SageJob"]
