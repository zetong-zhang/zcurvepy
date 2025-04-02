from setuptools import Extension, setup, find_packages

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
    url='https://zcurvepy-docs.readthedocs.io/',
    download_url='https://pypi.org/project/zcurvepy/',
    description='High performance Python toolkit for the Z-curve theory',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords=['Bioinformatics', 'Z-curve', 'Machine Learning'],
    license=license,
    version='1.5.10',
    install_requires=requirements,
    python_requires='>=3.7, <=3.11',
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_data={"ZCurvePy.assets.image": ["*.png", "*.jpg", "*.ico"]},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'zcurve-encoder=ZCurvePy.RunZCurveEncoder:main',
            'zcurve-plotter=ZCurvePy.RunZCurvePlotter:main',
            'zcurve-plotter-3d=ZCurvePy.RunZCurvePlotter3D:main',
            'zcurve-segmenter=ZCurvePy.RunZCurveSegmenter:main',
        ],
    },
    ext_modules=[Extension("_ZCurvePy", sources=[
        "src/ZCurvePy.cpp", 
        "src/ZCurvePyCore.cpp", 
        "src/ZCurvePyAPIs.cpp"
    ])]
)
