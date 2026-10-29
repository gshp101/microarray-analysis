"""Curated public GEO datasets used to develop and validate the pipeline.

Each entry says how to split the samples into a case group and a control
group. The regexes are matched against each sample's title, source and
characteristics text. `expected` holds the group sizes reported by the
original study, so a load can be checked against them; rules marked
verified=False have not yet been checked against the downloaded metadata.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class DatasetSpec:
    key: str
    accession: str
    description: str
    case: str  # name of the case group
    control: str  # name of the control group
    rules: dict[str, str]
    expected: dict[str, int]
    platform: str | None = None
    sample_filter: str | None = None
    exclude: str | None = None
    verified: bool = False
    reference: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)


DATASETS: dict[str, DatasetSpec] = {
    spec.key: spec
    for spec in [
        DatasetSpec(
            key="alzheimers_brain",
            accession="GSE5281",
            description="Alzheimer's disease vs non-demented control brain, laser-captured "
                        "neurons from six regions (Affymetrix HG-U133 Plus 2.0)",
            case="alzheimers", control="control",
            rules={"alzheimers": r"alzheimer|\baffected\b", "control": r"\bnormal\b|\bcontrol\b"},
            expected={"alzheimers": 87, "control": 74},
            platform="GPL570",
            verified=True,
            reference="Liang et al. 2008, Physiol Genomics; PMID 18270320",
            tags=("neurodegeneration", "dementia"),
        ),
        DatasetSpec(
            key="lung_adenocarcinoma",
            accession="GSE10072",
            description="Lung adenocarcinoma vs non-tumor lung tissue from smokers and "
                        "never smokers (Affymetrix HG-U133A)",
            case="tumor", control="normal",
            # titles look like "Lung Tumor_GT00006" / "Normal Lung_GT00006"
            rules={"tumor": r"adenocarcinoma", "normal": r"\bnormal\b"},
            expected={"tumor": 58, "normal": 49},
            platform="GPL96",
            verified=True,
            reference="Landi et al. 2008, PLoS One; PMID 18297132",
            tags=("cancer", "lung"),
        ),
        DatasetSpec(
            key="breast_cancer",
            accession="GSE15852",
            description="Breast tumor vs paired adjacent normal tissue "
                        "(Affymetrix HG-U133A)",
            case="tumor", control="normal",
            rules={"tumor": r"\bcancer\b|\btumou?r\b", "normal": r"\bnormal\b"},
            expected={"tumor": 43, "normal": 43},
            platform="GPL96",
            verified=True,
            reference="Ni et al. 2010, Breast Cancer Res Treat; PMID 20097481",
            tags=("cancer", "breast"),
        ),
        DatasetSpec(
            key="colorectal_cancer",
            accession="GSE32323",
            description="Colorectal cancer vs paired normal colon tissue; cell line "
                        "samples excluded (Affymetrix HG-U133 Plus 2.0)",
            case="cancer", control="normal",
            rules={"cancer": r"\bcancer\b|\btumou?r\b", "normal": r"\bnormal\b"},
            expected={"cancer": 17, "normal": 17},
            platform="GPL570",
            exclude=r"cell line",
            verified=True,
            reference="Khamas et al. 2012, Anticancer Res; PMID 23155247",
            tags=("cancer", "colorectal"),
        ),
    ]
}


def get(key_or_accession: str) -> DatasetSpec:
    if key_or_accession in DATASETS:
        return DATASETS[key_or_accession]
    for spec in DATASETS.values():
        if spec.accession.upper() == key_or_accession.upper():
            return spec
    raise KeyError(f"unknown dataset {key_or_accession!r}; known: {', '.join(DATASETS)}")
