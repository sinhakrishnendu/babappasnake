from pathlib import Path
import shutil

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.sdist import sdist as _sdist


ROOT = Path(__file__).resolve().parent
BUNDLE_DIR = ROOT / "src" / "babappasnake" / "_bundle"
BUNDLE_SOURCES = {
    "Snakefile": ROOT / "Snakefile",
    "workflow": ROOT / "workflow",
    "envs": ROOT / "envs",
    "config": ROOT / "config",
}


def refresh_bundle():
    BUNDLE_DIR.mkdir(parents=True, exist_ok=True)
    init_file = BUNDLE_DIR / "__init__.py"
    if not init_file.exists():
        init_file.write_text('"""Bundled Snakemake workflow assets."""\n')

    for path in BUNDLE_DIR.iterdir():
        if path.name == "__init__.py":
            continue
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()

    for name, source in BUNDLE_SOURCES.items():
        destination = BUNDLE_DIR / name
        if source.is_dir():
            ignore = shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo")
            if name == "config":
                ignore = shutil.ignore_patterns(
                    "__pycache__",
                    "*.pyc",
                    "*.pyo",
                    "package_live_override.yaml",
                )
            shutil.copytree(
                source,
                destination,
                ignore=ignore,
            )
        else:
            shutil.copy2(source, destination)


class build_py(_build_py):
    def run(self):
        refresh_bundle()
        super().run()


class sdist(_sdist):
    def run(self):
        refresh_bundle()
        super().run()


setup(cmdclass={"build_py": build_py, "sdist": sdist})
