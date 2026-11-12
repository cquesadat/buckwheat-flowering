from setuptools import find_packages, setup

setup(
    name='BuckwheatSeg',
    version='0.1.0',
    description='Mask R-CNN instance segmentation of buckwheat flowers',
    packages=find_packages(),
    zip_safe=False,
    python_requires='>=3.7,<3.8',
    install_requires=[
        'tensorflow-gpu==1.15.5',
        'tensorpack',
        'numpy<1.20',
        'opencv-python',
        'pycocotools',
        'scipy',
        'scikit-image',
        'Pillow',
        'Shapely',
        'tqdm',
        'tabulate',
        'termcolor',
        'six',
    ],
)
