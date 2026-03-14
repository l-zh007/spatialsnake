from setuptools import setup, find_packages, Command
import os

with open("requirements.txt") as req:
    requirements = [line.strip() for line in req if line.strip() and not line.startswith('#')]

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

with open("README.md", "r", encoding="utf-8") as readme:  # 添加encoding
    long_description = readme.read()

setup(
    name="spatialsnake",
    version="0.1.0",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.12.11",
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['spatialsnake=spatialsnake.command_line:main'],
    },
    author="lzh",
    author_email="1714074171@qq.com",
    description="spatialsnake",
    keywords="spatial transcription RNA analysis",
    cmdclass={'clean': CleanCommand},
    url="https://github.com/l-zh007/spatialsnake",
    classifiers=[
        'License :: OSI Approved :: MIT License',
        "Programming Language :: Python :: 3",
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
