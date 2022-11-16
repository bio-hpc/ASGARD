"""
This package supports the Amber updating script. It was written to automatically
check for, download, and/or apply bug fixes and updates for Amber.

This program was conceptualized and written by Jason Swails from the research
group of Adrian Roitberg at the University of Florida. It is released pursuant
to the GPL license, which can be found bundled inside AmberTools.

This is a rewrite of the original version, which was expanded and improved based
on the experiences I gained from the last version. My hope is to make it more
generalized and fix some of the flaws of the previous version.  I hope to better
support proxies, for instance, and to provide download progress for large
patches, etc.
"""

__author__ = 'Jason M. Swails'
__version__ = '14.0'
__license__ = 'GPL v3'

__all__ = ['downloader', 'exceptions', 'main', 'patch', 'patchlist',
           'preferences', 'progressbar', 'repos', 'utils', 'upgrade']
