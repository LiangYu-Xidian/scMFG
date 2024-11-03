from setuptools import setup
from setuptools import find_packages


def setup_package():
    metadata = dict(
        name='scmfg',
        version="1.0",
        description='scMFG: a single-cell Multi-omics Integration Method based on Feature Grouping',
        url='https://github.com/LiangYu-Xidian/scMFG',
        author='Yu',
        long_description=open('README.md').read(),
        long_description_content_type='text/markdown',  # 如果README是Markdown格式
        packages=find_packages(),
        install_requires=open('requirements.txt').read().splitlines()
    )
    setup(**metadata)


setup_package()
