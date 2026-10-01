"""Which GDC files a cohort's MAFs are drawn from.

A :class:`FileSpec` is the set of ``/files`` conditions that define one kind of
somatic mutation file. :data:`WXS_ENSEMBLE_MAF` is the open-access exome
ensemble MAF used by the documented runs; other file kinds are one constant
away, without touching the query code.
"""

from dataclasses import dataclass, field

from .api import gdc_and, gdc_in


@dataclass(frozen=True)
class FileSpec:
    """A set of GDC ``/files`` conditions identifying one kind of MAF.

    Every field other than ``name`` is a tuple of accepted values for the
    corresponding GDC field, or ``None`` to leave that field unconstrained.
    All constrained fields must hold for the *same* file, which is why
    :meth:`filter` is used against ``/files`` rather than ``/cases``.

    Parameters
    ----------
    name : str
        Short slug used in output file names, e.g. ``"wxs_ensemble_maf"``.
    title : str
        Human-readable name used in log messages and report text,
        e.g. ``"open WXS ensemble MAF"``. Defaults to ``name`` with
        underscores replaced by spaces.
    description : str
        One line describing the file kind, for logs and provenance records.
    data_category : tuple of str or None
        Accepted values of the GDC ``data_category`` field.
    data_type : tuple of str or None
        Accepted values of ``data_type``.
    experimental_strategy : tuple of str or None
        Accepted values of ``experimental_strategy``.
    workflow_type : tuple of str or None
        Accepted values of ``analysis.workflow_type``.
    access : tuple of str or None
        Accepted values of ``access``; ``("open",)`` keeps the run
        token-free.
    extra : tuple of (str, tuple) pairs
        Any further ``field, values`` conditions to require.
    """

    name: str
    title: str = ""
    description: str = ""
    data_category: tuple[str, ...] | None = None
    data_type: tuple[str, ...] | None = None
    experimental_strategy: tuple[str, ...] | None = None
    workflow_type: tuple[str, ...] | None = None
    access: tuple[str, ...] | None = None
    extra: tuple[tuple[str, tuple[str, ...]], ...] = field(default_factory=tuple)

    def __post_init__(self):
        """Default ``title`` to ``name`` with underscores replaced by spaces."""
        if not self.title:
            object.__setattr__(self, "title", self.name.replace("_", " "))

    def clauses(self):
        """Return the spec's filter clauses, without any case restriction."""
        fields = [
            ("data_category", self.data_category),
            ("data_type", self.data_type),
            ("experimental_strategy", self.experimental_strategy),
            ("analysis.workflow_type", self.workflow_type),
            ("access", self.access),
        ]
        clauses = [
            gdc_in(name, values) for name, values in fields if values is not None
        ]
        clauses += [gdc_in(name, values) for name, values in self.extra]
        return clauses

    def filter(self, case_ids):
        """Return the ``/files`` filter for this spec restricted to ``case_ids``."""
        return gdc_and(gdc_in("cases.case_id", case_ids), *self.clauses())

    def as_record(self):
        """Return the spec as a JSON-serialisable dict for provenance records."""
        return {
            "name": self.name,
            "title": self.title,
            "description": self.description,
            "data_category": self.data_category,
            "data_type": self.data_type,
            "experimental_strategy": self.experimental_strategy,
            "workflow_type": self.workflow_type,
            "access": self.access,
            "extra": {name: list(values) for name, values in self.extra},
        }


#: Open-access, GRCh38 ``Masked Somatic Mutation`` files from the Aliquot
#: Ensemble Somatic Variant Merging and Masking workflow: one MAF per
#: tumour-normal aliquot pair, exome (WXS). No token is needed to download
#: these. This is the spec the documented pipeline runs use.
WXS_ENSEMBLE_MAF = FileSpec(
    name="wxs_ensemble_maf",
    title="open WXS ensemble MAF",
    description="open-access WXS ensemble MAF (Masked Somatic Mutation, GRCh38)",
    data_category=("Simple Nucleotide Variation",),
    data_type=("Masked Somatic Mutation",),
    experimental_strategy=("WXS",),
    workflow_type=("Aliquot Ensemble Somatic Variant Merging and Masking",),
    access=("open",),
)
