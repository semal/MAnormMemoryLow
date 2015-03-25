#!/usr/bin/env python

import sys
from distutils.core import setup, Extension


def main():
    if not float(sys.version[:3]) >= 2.5:
        sys.stderr.write(
            "CRITICAL: Python version must be greater than or equal to 2.5! python 2.6.1 is recommended!\n")
        sys.exit(1)
    setup(name="MAnorm-Package",
          version="0.0.1",
          description="MAnorm, Memory low",
          author='semal',
          author_email='gongzhaohui@picb.ac.cn',
          url='http://www.github.com/semal/',
          package_dir={'MAnormMemoryLow': 'lib'},
          packages=['MAnormMemoryLow'],
          scripts=['bin/MAnormMemoryLow'],

          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Environment :: Web Environment',
              'Intended Audience :: Developers',
              'License :: OSI Approved :: Artistic License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: Microsoft :: Windows',
              'Operating System :: POSIX',
              'Programming Language :: Python',
              'Topic :: Database',
          ],)


if __name__ == '__main__':
    main()

