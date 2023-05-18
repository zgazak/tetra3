import os
from setuptools import setup, find_packages


def readme():
    with open("README.rst") as f:
        return f.read()


def read(rel_path):
    # type: (str) -> str
    here = os.path.abspath(os.path.dirname(__file__))
    # intentionally *not* adding an encoding option to open, See:
    #   https://github.com/pypa/virtualenv/issues/201#issuecomment-3145690
    with open(os.path.join(here, rel_path)) as fp:
        return fp.read()


def get_version(rel_path):
    # type: (str) -> str
    for line in read(rel_path).splitlines():
        if line.startswith("__version__"):
            # __version__ = "0.9"
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    raise RuntimeError("Unable to find version string.")


setup(
    name="sensor_capsule",
    version=get_version("sensor_capsule/__init__.py"),
    description="Collection Request and Response / Sensor interface",
    long_description=readme(),
    author="Zach Gazak",
    author_email="zachgazak@odysseyconsult.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "flask",
        "requests",
        "numpy",
        "pytz",
        "pathlib",
        "astropy",
    ],
    entry_points={
        "console_scripts": [
            "system=sensor_capsule.deploy.system.__main__:main",
            "bolt=sensor_capsule.deploy.bolt_app.__main__:main",
            "spout=sensor_capsule.deploy.spout_app.__main__:main",
            "registry=sensor_capsule.deploy.registry_app.__main__:main",
            "bolt_pipe=sensor_capsule.deploy.bolt_pipe.__main__:main",
            "spout_pipe=sensor_capsule.deploy.spout_pipe.__main__:main",
        ]
    },
    dependency_links=[],
    test_suite="nose.collector",
    tests_require=["nose"],
    zip_safe=False,
)
