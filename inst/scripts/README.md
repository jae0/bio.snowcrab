Most scripts now exist as markdown documented scripts in
[../markdown](../markdown). The numerical prefix indicates the sequence that
should be followed:

- 0X.*   : core assessment (required)

- 1X.*   : other research 

- 99*  : deprecated or testing scripts.


The file: [../markdown/Makefile](../markdown/Makefile) gives the recipes for creating html, latex, pdf, etc. for reports and summaries. The recipes us gnu-make which gives the commands required. VScode has a Makefile pluggin that can runthings for youif you have pandoc/quarto and gnu-make already installed. Alternatively, use vscode or r-studio to render them internally (but that requires you configure them appropriately -- this is a moving target so check current online resources). 
