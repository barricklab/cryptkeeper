from setuptools import setup, find_packages

# Gets Dependencies
with open('requirements.txt', 'r') as requirements_file:
    requirements = []
    for line in requirements_file:
        requirements.append(line.strip())

setup(
    name="cryptkeeper",
    version="0.1",
    packages={"": "bin"},

    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    install_requires=requirements,

    # metadata to display on PyPI
    author="Barrick Lab",
    author_email="jbarrick@cm.utexas.edu",
    description="Package for predicting functional components of DNA sequences",
    keywords="DNA, RNA, Cryptkeeper, Barrick, RBS, Biology, Bio",
    url="https://github.com/barricklab/cryptkeeper",   # project home page, if any
    project_urls={
        "Bug Tracker": "https://github.com/barricklab/cryptkeeper/issues",
        "Source Code": "https://github.com/barricklab/cryptkeeper",
    },
    classifiers=[
        "License :: GNU General Public License v2.0"
    ]

    # could also include long_description, download_url, etc.
)