```
Usage: ./phi_analysis.R [options] file1 file2 ...


Options:
    --phi=CHARACTER
    comma separated string defining phi values: [default= 0.177,0.0442]

    -c CHARACTER, --col=CHARACTER
	collection table with name and collection columns

    -h, --help
    Show this help message and exit
```

```
Usage: ./duplicate_analisys.R [options] file1 file2 ...

Options:
	-t OPTION, --phi=OPTION
		phi threshold to infer duplication [default= 0.45]

	-v OPTION, --vertexinfo=OPTION
		tab delimited file, minimum 2 columns: name collection

	-r OPTION, --resultsdir=OPTION
		results directory [default= results]

	-h, --help
		Show this help message and exit
```

```
Usage: ./fastivs_analysis.R [options] file


Options:
	-t OPTION, --phi=OPTION
		phi threshold to infer independence [default= 0.0442]

	-n OPTION, --nreq=OPTION
		number of ivs requested to fastindep [default= 2]

	-m OPTION, --nread=OPTION
		number of ivs to read from fastindep [default= 2]

	-v OPTION, --vertexinfo=OPTION
		tab delimited file, minimum 2 columns: name collection

	-r OPTION, --resultsdir=OPTION
		results directory [default= fastivs]

	-h, --help
		Show this help message and exit
```
