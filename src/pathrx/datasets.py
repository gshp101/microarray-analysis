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
    lineage: str | None = None  # tissue of matching cancer cell lines in DepMap/CCLE, e.g. "LUNG"


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
            lineage="LUNG",
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
            lineage="BREAST",
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
            lineage="LARGE INTESTINE",
        ),
        DatasetSpec(
            key="alzheimers_aging_brain",
            accession="GSE48350",
            description="Alzheimer's disease vs cognitively intact controls aged 60 and over, "
                        "four brain regions (hippocampus, entorhinal cortex, superior frontal "
                        "and postcentral gyrus); controls under 60 left out so age does not "
                        "drive the contrast (Affymetrix HG-U133 Plus 2.0)",
            case="alzheimers", control="control",
            # AD titles look like "hippocampus_female_85_AD_12", controls like
            # "Hippocampus_male_45yrs_indiv12" (the series design labels AD cases "AD")
            rules={"alzheimers": r"_AD_", "control": r"_indiv\d+"},
            # study: 80 AD and 173 controls aged 20-99, of whom 93 are 60 or older
            expected={"alzheimers": 80, "control": 93},
            platform="GPL570",
            sample_filter=r"age \(yrs\): (?:[6-9]\d|1\d\d)\b",
            verified=True,
            reference="Berchtold et al. 2008, PNAS; PMID 18832152 (AD cases also in Cribbs et al. "
                      "2012, J Neuroinflammation; PMID 22824372)",
            tags=("neurodegeneration", "dementia", "aging"),
        ),
        DatasetSpec(
            key="liver_cancer",
            accession="GSE14520",
            description="Hepatocellular carcinoma vs paired non-tumor liver tissue, main "
                        "cohort (Affymetrix HT HG-U133A, GPL3921)",
            case="tumor", control="nontumor",
            # titles "Liver Tumor Tissue LCS-079A" / "Liver Non-Tumor Tissue LCS-079B"
            rules={"tumor": r"liver tumor tissue", "nontumor": r"non-tumor tissue"},
            expected={"tumor": 225, "nontumor": 220},
            platform="GPL3921",
            verified=True,
            reference="Roessler et al. 2010, Cancer Res; PMID 21159642",
            tags=("cancer", "liver"),
            lineage="LIVER",
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
