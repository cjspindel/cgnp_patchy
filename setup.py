from setuptools import setup


setup(
    # Self-descriptive entries which should always be present
    name='cgnp_patchy',
    author='Caroline Spindel',
    author_email='cjs323@lehigh.edu',
    license='MIT',
    version='0.0.0',
    description='An mBuild recipe for generating parameterized models of polymer-tethered, coarse-grained silica nanoparticles.',
    zip_safe=False,
    entry_points={
        'mbuild.plugins':[
        "cgnp_patchy = cgnp_patchy.cgnp_patchy:cgnp_patchy"
        ]
        }
    )
