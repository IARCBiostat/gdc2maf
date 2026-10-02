"""Run records: where the output goes, and what produced it.

A library should not grab hold of the root logger, so nothing here happens on
import. Call :func:`configure_logging` to see the package's progress messages
(the CLI does this for you), and :func:`log_environment` to record the
interpreter and package versions that produced a run.

:func:`tee_stdout` is the heavier option: it redirects ``sys.stdout`` and
``sys.stderr`` so that output from *other* libraries — a progress bar in a
downstream tool, say — also lands in the log file.
"""

import logging
import platform
import sys
from datetime import datetime
from importlib import metadata

logger = logging.getLogger(__name__)

#: Distributions whose versions are recorded by :func:`log_environment`.
DEFAULT_PACKAGES = ["gdc2maf", "pandas", "numpy", "requests"]

PLAIN_FORMAT = "%(message)s"

#: Sentinel for "whatever ``sys.stdout`` is when the call is made". A plain
#: ``stream=sys.stdout`` default would capture the stdout in place at import
#: time and so miss a later :func:`tee_stdout`.
USE_STDOUT = object()


def configure_logging(
    log_path=None, level=logging.INFO, stream=USE_STDOUT, fmt=PLAIN_FORMAT, mode="w"
):
    """Send the package's log messages to ``stream`` and optionally to a file.

    Attaches handlers to the ``gdc2maf`` logger only, so an application's own
    logging configuration is left alone. Calling it again replaces the
    handlers it added before, which makes it safe to call from a notebook.

    Parameters
    ----------
    log_path : str or None, optional
        File to write the log to. ``None`` logs to ``stream`` only.
    level : int, optional
        Logging level for the package's messages.
    stream : file object, None or USE_STDOUT, optional
        Stream to log to. The default resolves to ``sys.stdout`` as it is at
        the time of the call, so logging follows a :func:`tee_stdout` installed
        beforehand. ``None`` logs to the file only.
    fmt : str, optional
        Logging format. The default prints the message alone, so the log reads
        like the pipeline's own output.
    mode : str, optional
        File mode; ``"w"`` (the default) overwrites the log on every run,
        ``"a"`` appends.

    Returns
    -------
    logging.Logger
        The configured ``gdc2maf`` logger.
    """
    package_logger = logging.getLogger("gdc2maf")
    for handler in [
        h for h in package_logger.handlers if getattr(h, "_gdc2maf", False)
    ]:
        package_logger.removeHandler(handler)
        handler.close()

    if stream is USE_STDOUT:
        stream = sys.stdout

    formatter = logging.Formatter(fmt)
    handlers = []
    if stream is not None:
        handlers.append(logging.StreamHandler(stream))
    if log_path is not None:
        handlers.append(logging.FileHandler(log_path, mode=mode))
    for handler in handlers:
        handler.setFormatter(formatter)
        handler._gdc2maf = True
        package_logger.addHandler(handler)

    package_logger.setLevel(level)
    package_logger.propagate = False
    if log_path is not None:
        logger.info(
            f"Run started {datetime.now().isoformat(timespec='seconds')}: "
            f"{' '.join(sys.argv)}"
        )
    return package_logger


class _Tee:
    """File-like object writing to several streams.

    The first stream is the "real" one: queries about the underlying
    terminal (``fileno``, ``isatty``, ``encoding``) are answered by it.
    Libraries such as ``alive_progress`` (used by SigProfilerExtractor)
    inspect the stream that way and fail on objects that only write.
    """

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
        return len(text)

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        for stream in self.streams:
            stream.flush()

    def fileno(self):
        return self.streams[0].fileno()

    def isatty(self):
        return self.streams[0].isatty()

    @property
    def encoding(self):
        return getattr(self.streams[0], "encoding", "utf-8")

    @property
    def errors(self):
        return getattr(self.streams[0], "errors", None)

    def writable(self):
        return True

    def readable(self):
        return False

    def seekable(self):
        return False

    @property
    def closed(self):
        return False


def tee_stdout(log_path, mode="w"):
    """Copy ``sys.stdout`` and ``sys.stderr`` to ``log_path`` as well.

    Use this when output from libraries outside this package has to reach the
    log file too; for the package's own messages, :func:`configure_logging` is
    enough. Errors (tracebacks) are included because stderr is copied too.

    Parameters
    ----------
    log_path : str
        Path of the log file.
    mode : str, optional
        File mode; ``"w"`` (the default) overwrites it on every run.

    Returns
    -------
    file object
        The open log file; closed at interpreter exit.
    """
    log_file = open(log_path, mode, buffering=1)
    sys.stdout = _Tee(sys.__stdout__, log_file)
    sys.stderr = _Tee(sys.__stderr__, log_file)
    return log_file


def environment_info(packages=None):
    """Return the interpreter, platform and package versions as a dict.

    Parameters
    ----------
    packages : list of str or None, optional
        Distribution names to report; missing ones are reported as
        ``not installed``. ``None`` uses :data:`DEFAULT_PACKAGES`.

    Returns
    -------
    dict
        ``python``, ``executable``, ``platform`` and one entry per package.
    """
    versions = {}
    for package in packages if packages is not None else DEFAULT_PACKAGES:
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return {
        "python": platform.python_version(),
        "executable": sys.executable,
        "platform": platform.platform(),
        **versions,
    }


def log_environment(packages=None):
    """Log the interpreter, platform and package versions of this run.

    Parameters
    ----------
    packages : list of str or None, optional
        Distribution names to report; ``None`` uses :data:`DEFAULT_PACKAGES`.
        Pass extra names to record a downstream tool's version in the same
        place, e.g. ``DEFAULT_PACKAGES + ["SigProfilerExtractor"]``.

    Returns
    -------
    dict
        The dict returned by :func:`environment_info`.
    """
    info = environment_info(packages)
    logger.info("Environment:")
    for key, value in info.items():
        logger.info(f"  {key}: {value}")
    return info
