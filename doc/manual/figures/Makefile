inkscape = inkscape --export-area-drawing

random-walk.png: random-walk.svg
	$(inkscape) --export-png=$@ --export-width 500 $^

random-walk.pdf: random-walk.svg
	$(inkscape) --export-pdf=$@ $^
