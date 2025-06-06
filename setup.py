from setuptools import Extension, setup, find_packages
import numpy as np

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("LICENSE", "r", encoding="utf-8") as fh:
    license = fh.read()

with open("requirements.txt", encoding="utf-8") as fh:
    requirements = fh.readlines()

setup(
    name='zcurvepy',
    author='TUBIC',
    author_email='fgao@tju.edu.cn',
    maintainer='Zetong Zhang',
    maintainer_email='zhangzetong@tju.edu.cn',
    url='https://zcurvehub-docs.readthedocs.io/',
    download_url='https://pypi.org/project/zcurvepy/',
    description='High performance Python toolkit for the Z-curve theory',
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    keywords=['Bioinformatics', 'Z-curve', 'Machine Learning'],
    license=license,
    version='1.6.0',
    install_requires=requirements,
    python_requires='>=3.7, <=3.11',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={"ZcurvePy.assets.image": ["*.png", "*.jpg", "*.ico"]},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'zcurve-encoder=ZcurvePy.Scripts.RunZcurveEncoder:main',
            'zcurve-plotter=ZcurvePy.Scripts.RunZcurvePlotter:main',
            'zcurve-plotter-3d=ZcurvePy.Scripts.RunZcurvePlotter3D:main',
            'zcurve-segmenter=ZcurvePy.Scripts.RunZcurveSegmenter:main',
        ],
    },
    ext_modules=[Extension("_ZcurvePy", 
        include_dirs=[np.get_include()],
        sources=[
            "ZcurvePy/include/ZcurvePy.cpp", 
            "ZcurvePy/include/ZcurvePyCore.cpp", 
            "ZcurvePy/include/ZcurvePyAPIs.cpp",
            "ZcurvePy/include/ZislandFinder.cpp"
        ]
    )]
)
