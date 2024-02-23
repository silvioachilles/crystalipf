from distutils.core import setup

setup(
    name='crystalipf',
    version='0.0.1',
    packages=['crystalipf'],
    author='Silvio Achilles',
    author_email='achilles.develop@proton.me',
    description='Easy-to-use library for calculating the rgb values of the IPF map of cubic materials.',
    install_requires=["numpy", "matplotlib"]
)
