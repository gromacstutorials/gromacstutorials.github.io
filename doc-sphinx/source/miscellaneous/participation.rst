.. contribute-label:

Wanna contribute?
*****************

Your feedback is always useful. You can also suggest idea for new a tutorial, or a new section to an existing tutorial, just keep in mind that:

- the content must be of interest for everyone, not just for the specific issue you are trying to solve,
- the tutorials must remain relatively simple. 

As a tester
===========

Report broken link and typo by `email`_ or by posting a new issue on |github-issue|.

.. _email: simon.gravelle@live.fr

.. |github-issue| raw:: html

   <a href="https://github.com/gromacstutorials/gromacstutorials.github.io/issues" target="_blank">Github</a>


As a writer
===========

Propose a new tutorial or a modification to an existing tutorial.
To do so, fork the gromacstutorials repository on github, make your changes,
and submit a |pull_request|.

.. |pull_request| raw:: html

   <a href="https://github.com/gromacstutorials/gromacstutorials.github.io/pulls" target="_blank">pull request</a>

Build *gromacstutorial* locally
===============================

*gromacstutorial* can be build locally on your computer using sphinx. Fork
the repository by typing in a terminal:

..  code-block:: bash

    git clone https://github.com/gromacstutorials/gromacstutorials.github.io.git

Then, go to docs/sphinx/, and execute:

..  code-block:: bash

    bash build.sh

The docs/index.html can be opened with a web browser.
The tutorials are located in docs/sphinx/source/tutorials and written in rst format. 
For each tutorials, the corresponding gromacs input files are located here: docs/inputs.
