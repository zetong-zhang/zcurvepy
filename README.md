# ZcurvePy (Lastest version 1.6.0)

*A high performance Python toolkit for the Z-curve theory developed by **T**ianjin **U**niversity **B**io**I**nformatics **C**enter (**TUBIC**)*  

## Note:  
This page only provides a brief description of the software.  
For more information, please visit our full API and help documentation on [this website](https://zcurvehub-docs.readthedocs.io).

## Contents
- [Overview](#title1)
    - [Core Capabilities](#subtitle1)
    - [Technical Highlights](#subtitle2)
- [Installation](#title2)
    - [Python Requirements](#subtitle3)
    - [Operating System](#subtitle4)
- [Quickstart](#title3)
    - [Generate Z-curves and derivatives](#section1)
    - [Extract Z-curve parameter features](#section2)
    - [Segmentation for DNA sequences](#section3)
    - [Build a simple gene recognizer](#section4)
- [Web Server](#title4)
- [Reference](#title5)
- [Citing](#title6)
- [Contact](#title7)
- [Acknowledgement](#title8)

## Overview <a id="title1"></a>
ZcurvePy aims to build a comprehensive, one-stop bioinformatics platform designed to streamline nucleic acid sequence analysis through the mathematical framework of the Z-curve theory. It empowers researchers to extract DNA/RNA structural features, identify functional genomic elements, and build predictive models with cutting-edge computational tools.
### Core Capabilities <a id="subtitle1"></a>
- **Genome Visualization with Z-curves**  
Generate Z-curves and their derivatives (e.g., GC disparity, CpG profile) from raw sequences, with customizable visualization of geometric features displayed in 2D or 3D. Supports FASTA and GenBank formats. 
- **Feature Extraction and Selection**  
Extract and select features using Z-curve parameters more customarily and flexibly compared to non-standalone modules integrated into other software, and explores its powerful application in gene prediction, promoter classification, replication origin recognition, etc. with machine learning or deep learning.
- **Accurate Curve Segmentation**  
Detect critical structural boundaries using genome order index algorithm, identifying candidate regions for replication origins, horizontal gene transfer event or CpG islands in eukaryotic genomes.
- **Build Classification Models**  
Construct nucleic acid sequence classifier with biological function based on machine learning or deep learning framework, high-precision protein gene recognizers for specific species taxa of prokaryotes, which is very useful when studying newly sequenced or resequenced species that are closely related.
### Technical Highlights <a id="subtitle2"></a>
1. **High-Performance Hybrid Architecture**
    - **C/C++ Acceleration**  
    Core algorithmic modules are implemented natively in C/C++ and seamlessly integrated with Python via dynamic libraries (DLL/SO), where C++ classes and functions are wrapped into Python-callable objects using native Python C/C++ APIs, balancing development efficiency with runtime performance.
    - **Parallel Computing**  
    Allows multi-threaded parallelization, achieving 4-6x speedup for large-scale genomic data processing (e.g., 765-bit Z-curve parameters generation for *S. cerevisiae*'s CDS sequences takes 0.3 seconds vs. 1.3 seconds in single-threaded mode)
2. **Cross-Paradigm Interfaces**
    - **Command-Line Interface**  
    Streamlined CLI commands for batch processing and pipeline integration, ideal for bioinformatics workflows, e.g.
        ```bash
        zcurve-encoder -f example.fa -s settings.json -o features.csv
        ```
    - **Python API**  
    Object-oriented interfaces for developers, enabling customizable workflows and real-time result callbacks, e.g.
        ```python
        # Init ZcurveEncoder
        from ZcurvePy import BatchZcurveEncoder
        hyper_params = [ ... ]
        encoder = BatchZcurveEncoder(hyper_params, n_jobs=8)
        # Load and process data
        from Bio.SeqIO import parse
        records = parse("example.fa", "fasta")
        features = encoder(records)
        ```
3. **Ecosystem Integration**  
    - **Data Connectivity**  
    Built-in integration with [Biopython](https://pypi.org/project/biopython/) and [Entrez](https://www.ncbi.nlm.nih.gov/sites/batchentrez) modules for direct sequence retrieval from NCBI databases (e.g., `download_acc("NC_000913")`), with automated parsing of FASTA/GenBank formats.
    - **ML Compatibility**  
    Extracted Z-curve features are directly compatible with [scikit-learn](https://scikit-learn.org/) (traditional ML) and [PyTorch](https://scikit-learn.org/) (deep learning), including pre-trained models (e.g., Ori-FinderH, Nmix).
    - **Visualization Tools**  
    Export Z-curve trajectories as [Matplotlib](https://matplotlib.org/) static plots or [Plotly](https://plotly.com/) interactive HTML.
## Installation <a id="title2"></a>
Python includes the package management system `pip` which should allow you to install ZcurvePy and its dependencies if needed, as well as upgrade or uninstall with just one command in the terminal:
```bash
python -m pip install zcurvepy
python -m pip install --upgrade zcurvepy

python -m pip uninstall zcurvepy
```
Starting from 1.5.11, the return value types of some frequently used API functions have been modified from Python list to Numpy ndarray. Therefore, please install Numpy 1.x before compiling and installing.
```bash
python -m pip install "numpy>=1.21,<2.0"
```
### Python Requirements <a id="subtitle3"></a>
Python 3.7, 3.8, 3.9, 3.10, 3.11, 3.12 are supported. We currently recommend using Python 3.9.6 (https://www.python.org/downloads/release/python-396/)
### Operating System <a id="subtitle4"></a>
Windows 10/11, macOS, and Linux running on x86_64 arch are supported. Note that some features of Matplotlib may not work on systems without a graphic interface (e.g., RedHat).
## Quickstart <a id="title3"></a>  
1. Generate Z-curves and derivatives<a id="section1"></a>  
   (1) Python API implementation:  

    ```python
    from ZcurvePy import ZcurvePlotter
    from ZcurvePy.Util import download_acc
    from Bio.SeqIO import read
    import matplotlib.pyplot as plt
    # Download genomes (Please wait for several seconds)
    save_paths = download_acc("NC_000854.2,NC_000868.1")
    ax3d = plt.figure().add_subplot(projection='3d')
    ax2d = plt.figure().add_subplot()
    for save_path in save_paths:
        record = read(save_path, "fasta")
        # Calculate components of smoothed Z-curve
        plotter = ZcurvePlotter(record)
        n, x, y, _ = plotter.z_curve(window=1000)
        zp, _ = plotter.z_prime_curve(window=1000, return_n=False)
        # Matplotlib 2D display
        ax2d.plot(n[::10], x[::10], label=f"{record.id} RY-disparity")
        ax2d.plot(n[::10], y[::10], label=f"{record.id} MK-disparity")
        ax2d.plot(n[::10], zp[::10], label=f"{record.id} Z'n curve")
        ax2d.legend()
        # Matplotlib 3D display
        ax3d.plot(x[::10], y[::10], zp[::10], label=f"{record.id} Z-curve")
        ax3d.legend()
    plt.show()
    ```
    (2) Commandline usage:  
    For plotting curves in 2D mode:
    ```bash
    zcurve-plotter -a NC_000854.2,NC_000868.1 -s settings.json -p curves.png -o curves_2d.json
    ```
    and for plotting curves in 3D mode:
    ```bash
    zcurve-plotter-3d -a NC_000854.2,NC_000868.1 -s settings.json -p curves.png -o curves_3d.json
    ```
    where the settings should be a JSON like:
    ```json
    {
        "plotter": [
            {
                "window": 10000, 
                "intv": 100,
                "curve2d": "RY,MK,ZP",
                "curve3d": "RY:MK:ZP"
            },
            {
                "window": 10000, 
                "intv": 100,
                "curve2d": "RY,MK,ZP",
                "curve3d": "RY:MK:ZP"
            }
        ]
    }
    ```
    Note that the Z-curve Plotter's 3D mode allows you to select 3 of the 11 curves as components of x,y, and z.  
    For more help information about the setting file, please enter:

    ```bash
    zcurve-plotter --help
    ```
2. Extract Z-curve parameter features<a id="section2"></a>  
    (1) Python API implementation:  

    ```python
    from ZcurvePy.Util import download_acc, extract_CDS
    import numpy as np
    save_paths = download_acc("NC_000854.2", __name__, "gb")
    records = extract_CDS(save_paths[0])

    from ZcurvePy import ZcurveEncoder
    def encoding(record):
        """ simple batch processing """
        encoder = ZcurveEncoder(record)
        # Calculate and concatenate 765-bit Z-curve transform
        feature = encoder.mononucl_phase_transform(freq=True)
        feature = np.concatenate((feature, encoder.dinucl_phase_transform(freq=True)))
        feature = np.concatenate((feature, encoder.trinucl_phase_transform(freq=True)))
        feature = np.concatenate((feature, encoder.k_nucl_phase_transform(k=4, freq=True)))
        return feature

    features = np.array([encoding(record) for record in records])
    print(features.shape)
    ```

    or use another more powerful API to implement multi-threading:

    ```python
    from ZcurvePy import BatchZcurveEncoder
    # Define the hyper-paramsfor 765-bit Z-curve transform
    hyper_params = [
        {"k": 1, "freq": True}  # Same as mononucl_phase_transform(freq=True)
        {"k": 2, "freq": True}  # Same as dinucl_phase_transform(freq=True)
        {"k": 3, "freq": True}  # Same as trinucl_phase_transform(freq=True)
        {"k": 4, "freq": True}  # Same as k_nucl_phase_transform(k=4, freq=True)
    ]
    encoder = BatchZcurveEncoder(hyper_params, n_jobs=8)
    features = encoder(records)
    ```
    (2) Commandline usage:
    ```bash
    zcurve-encoder -a NC_000854.2 -s settings.json -e True -o features.csv
    ```
    where the setting file should be a JSON like:
    ```json
    {
        "encoder": {
            "hyper_params": [
                {"k": 1, "freq": true},
                {"k": 2, "freq": true},
                {"k": 3, "freq": true},
                {"k": 4, "freq": true}
            ],
            "n_jobs": 8
        }
    }
    ```
    For more help information about the setting file, please enter:
    ```bash
    zcurve-encoder --help
    ```
3. Segmentation for DNA sequences<a id="section3"></a>  
    Python API implementation:
    ```python
    from ZcurvePy.Util import download_acc
    from Bio.SeqIO import read
    from ZcurvePy import ZcurveSegmenter, ZcurvePlotter
    import matplotlib.pyplot as plt
    # Download data
    path = download_acc("CP001956.1")
    record = read(path[0], "fasta")
    # Segmentation
    segmenter = ZcurveSegmenter(mode='WS', min_len=50000)
    seg_points = segmenter.run(record)
    # Calculate z' curve for visualization
    plotter = ZcurvePlotter(record)
    n, zp, _ = plotter.z_prime_curve()
    # Visualization
    for point, _ in seg_points:
        plt.axvline(point, color='red')
    plt.plot(n, zp)
    plt.show()
    ```
    Commandline usage:
    ```bash
    zcurve-segmenter -a CP001956.1 -m WS -l 50000 -o seg_points.csv -v True
    ```
4. Build a simple gene recognizer <a id="section4"></a>
    This is an example of training a gene recognition model for *Escherichia coli* :
    ```python
    # We recommend turning on the acceleration system on Intel platforms
    from sklearnex import patch_sklearn
    patch_sklearn()
    from ZcurvePy.Util import download_acc, extract_CDS
    from ZcurvePy import ZcurveBuilder
    from Bio.SeqIO import read, parse
    # Load positive dataset
    path = download_acc("NC_000913.3", __name__, "gb")[0]
    pos_dataset = extract_CDS(path)
    builder = ZcurveBuilder(standard=True, n_jobs=8)
    builder.fit(pos_dataset)
    # Some sample sequences
    records = parse("samples.fa", "fasta")
    results = builder.predict(records)
    ```
## Web Server & Database <a id="title4"></a>
A free, flexible and interactive **ZcurveHub** Web Service as well as updated **ZcurveDB** is available at http://tubic.tju.edu.cn/zcurve/.

## Reference <a id="title5"></a>
[1] &nbsp; Guo FB, Ou HY, Zhang CT. ZCURVE: a new system for recognizing protein-coding genes in bacterial and archaeal genomes. Nucleic Acids Res. 2003 Mar 15;31(6):1780-9. doi: 10.1093/nar/gkg254. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/12626720/)  

[2] &nbsp; Zhang CT, Zhang R. A nucleotide composition constraint of genome sequences. Comput Biol Chem. 2004 Apr;28(2):149-53. doi: 10.1016/j.compbiolchem.2004.02.002. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/15130543/)  

[3] &nbsp; Zhang CT, Gao F, Zhang R. Segmentation algorithm for DNA sequences. Phys Rev E Stat Nonlin Soft Matter Phys. 2005 Oct;72(4 Pt 1):041917. doi: 10.1103/PhysRevE.72.041917. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/16383430/)  

[4] &nbsp; Gao F, Zhang CT. GC-Profile: a web-based tool for visualizing and analyzing the variation of GC content in genomic sequences. Nucleic Acids Res. 2006 Jul 1;34(Web Server issue):W686-91. doi: 10.1093/nar/gkl040. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/16845098/)   

[5] &nbsp; Zhang R, Zhang CT. A Brief Review: The Z-curve Theory and its Application in Genome Analysis. Curr Genomics. 2014 Apr;15(2):78-94. doi: 10.2174/1389202915999140328162433. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/24822026/)  

[6] &nbsp; Hua ZG, Lin Y, Yuan YZ, Yang DC, Wei W, Guo FB. ZCURVE 3.0: identify prokaryotic genes with higher accuracy as well as automatically and accurately select essential genes. Nucleic Acids Res. 2015 Jul 1;43(W1):W85-90. doi: 10.1093/nar/gkv491. Epub 2015 May 14. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/25977299/)  

[7] &nbsp; Wang D, Lai FL, Gao F. Ori-Finder 3: a web server for genome-wide prediction of replication origins in Saccharomyces cerevisiae. Brief Bioinform. 2021 May 20;22(3):bbaa182. doi: 10.1093/bib/bbaa182. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/34020544/)  

[8] &nbsp; Lai FL, Gao F. GC-Profile 2.0: an extended web server for the prediction and visualization of CpG islands. Bioinformatics. 2022 Mar 4;38(6):1738-1740. doi: 10.1093/bioinformatics/btab864.  [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/34954794/)  

[9] Yin ZN, Lai FL, Gao F. Unveiling human origins of replication using deep learning: accurate prediction and comprehensive analysis. Brief Bioinform. 2023 Nov 22;25(1):bbad432. doi: 10.1093/bib/bbad432.  [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/38008420/)  

[10] Geng YQ, Lai FL, Luo H, Gao F. Nmix: a hybrid deep learning model for precise prediction of 2'-O-methylation sites based on multi-feature fusion and ensemble learning. Brief Bioinform. 2024 Sep 23;25(6):bbae601. doi: 10.1093/bib/bbae601. [[Pubmed]](https://pubmed.ncbi.nlm.nih.gov/39550226/) 

## Citation <a id="title6"></a>
The paper on this work has not yet been published. If you would like to cite this software in your work, please contact us to discuss alternatives.

Z Zhang, Y Lin, H Luo, F Gao. ZcurveHub: an updated large-scale Z curve knowledgebase with Scalable Genome Analysis Framework.

## Contact <a id="title7"></a>
The offical website of TUBIC: https://tubic.org/ | https://tubic.tju.edu.cn .  
If you have any questions about this software, please contact fgao@tju.edu.cn .  

Copyright © Tianjin University BioInformatics Center  
No. 92 Weijin Road Nankai District  
Tianjin, China, 300072  
Telephone: +86-22-27402697  

## Acknowledgement <a id="title8"></a>
The authors wishes to thank Chun-Ting Zhang, Academician of the Chinese Academy of Sciences, who proposed the Z-curve theory and his collaborators Ren Zhang, Lin-Lin Chen, Hong-Yu Ou and Feng-Biao Guo for their significant contributions to the development of the theory.  <p align="right">— Staff of TUBIC, 2025-06-02</p>