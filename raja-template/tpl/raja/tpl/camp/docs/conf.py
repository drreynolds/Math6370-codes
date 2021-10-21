# -*- coding: utf-8 -*-
#
# NanoGUI documentation build configuration file, created by
# sphinx-quickstart on Mon Aug 22 20:05:54 2016.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys
import os
import shlex
import textwrap

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('.'))

# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# [[[ begin extensions marker ]]]
# Tell Sphinx to use both the `breathe` and `exhale` extensions
extensions = [
    'breathe',
    'exhale'
]

# Setup the `breathe` extension
breathe_projects = { "camp": "./doxyoutput/xml" }
breathe_default_project = "camp"

# Setup the `exhale` extension
import textwrap
exhale_args = {
    ############################################################################
    # These arguments are required.                                            #
    ############################################################################
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "Library API",
    "doxygenStripFromPath":  "../include",
    ############################################################################
    # Suggested optional arguments.                                            #
    ############################################################################
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin": textwrap.dedent('''
        INPUT       = ../include
        EXCLUDE_SYMBOLS = detail detail::* internal internal::* test test::*
        CLANG_ASSISTED_PARSING = YES
        CLANG_OPTIONS = -std=c++14 -stdlib=libc++ -DCAMP_DOX
        PREDEFINED += CAMP_HOST_DEVICE=""
        PREDEFINED += CAMP_DEVICE=""
        PREDEFINED += CAMP_CONSTEXPR14=constexpr
        PREDEFINED += CAMP_DOX
    '''),
    ############################################################################
    # HTML Theme specific configurations.                                      #
    ############################################################################
    # Fix broken Sphinx RTD Theme 'Edit on GitHub' links
    # Search for 'Edit on GitHub' on the FAQ:
    #     http://exhale.readthedocs.io/en/latest/faq.html
    "pageLevelConfigMeta": ":github_url: https://github.com/llnl/camp",
    ############################################################################
    # Main library page layout example configuration.                          #
    ############################################################################
    # "afterTitleDescription": textwrap.dedent(u'''
    #     Welcome to the developer reference to Exhale Companion.  The code being
    #     documented here is largely meaningless and was only created to test
    #     various corner cases e.g. nested namespaces and the like.
    #
    #     .. note::
    #
    #         The text you are currently reading was fed to ``exhale_args`` using
    #         the :py:data:`~exhale.configs.afterTitleDescription` key.  Full
    #         reStructuredText syntax can be used.
    #
    #     .. tip::
    #
    #        Sphinx / Exhale support unicode!  You're ``conf.py`` already has
    #        it's encoding declared as ``# -*- coding: utf-8 -*-`` **by
    #        default**.  If you want to pass Unicode strings into Exhale, simply
    #        prefix them with a ``u`` e.g. ``u"👽😱💥"`` (of course you would
    #        actually do this because you are writing with åçćëñtß or
    #        non-English 寫作 😉).
    # '''),
    # "afterHierarchyDescription": textwrap.dedent('''
    #     Below the hierarchies comes the full API listing.
    #
    #     1. The text you are currently reading is provided by
    #        :py:data:`~exhale.configs.afterHierarchyDescription`.
    #     2. The Title of the next section *just below this* normally defaults to
    #        ``Full API``, but the title was changed by providing an argument to
    #        :py:data:`~exhale.configs.fullApiSubSectionTitle`.
    #     3. You can control the number of bullet points for each linked item on
    #        the remainder of the page using
    #        :py:data:`~exhale.configs.fullToctreeMaxDepth`.
    # '''),
    # "fullApiSubSectionTitle": "Custom Full API SubSection Title",
    # "afterBodySummary": textwrap.dedent('''
    #     You read all the way to the bottom?!  This text is specified by giving
    #     an argument to :py:data:`~exhale.configs.afterBodySummary`.  As the docs
    #     state, this summary gets put in after a **lot** of information.  It's
    #     available for you to use if you want it, but from a design perspective
    #     it's rather unlikely any of your users will even see this text.
    # '''),
    ############################################################################
    # Individual page layout example configuration.                            #
    ############################################################################
    # Example of adding contents directives on custom kinds with custom title
    "contentsTitle": "Page Contents",
    "kindsWithContentsDirectives": ["class", "file", "namespace", "struct"],
    # This is a testing site which is why I'm adding this
    "includeTemplateParamOrderList": True,
    ############################################################################
    # useful to see ;)
    "verboseBuild": True,
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'
# [[[ end extensions marker ]]]

# on_rtd is whether we are on readthedocs.org, this line of code grabbed from docs.readthedocs.org
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:  # only import and set the theme if we're building docs locally
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
# [[[ end theme marker ]]]

# NOTE: this is only here for me to auto-place sections of conf.py in the docs
#       but isn't needed in producion releases
html_theme = 'sphinx_rtd_theme'

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'CAMP'
copyright = u'2019, Tom Scogland'
author = u'Tom Scogland'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']

# The reST default role (used for this markup: `text`) to use for all
# documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []

# If true, keep warnings as "system message" paragraphs in the built documents.
#keep_warnings = False

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False
