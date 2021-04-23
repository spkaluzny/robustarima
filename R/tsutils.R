"tssub" <- 
##
## given start and end, return td[start:end]
## 
function(td, start = NULL, end = NULL, ...)
{
	idx <- NULL
	if(!is.null(start)) {
		if(is(start, "timeDate")) {
			idx <- (td >= start)
		}
		else {
			idx <- (td >= timeDate(start, ...))
		}
	}
	if(!is.null(end)) {
		if(is(end, "timeDate")) {
			if(!is.null(idx)) {
				idx <- idx & (td <= end)
			}
			else {
				idx <- (td <= end)
			}
		}
		else {
			if(!is.null(idx)) {
				idx <- idx & (td <= timeDate(end, ...))
			}
			else {
				idx <- (td <= timeDate(end, ...))
			}
		}
	}
	if(!is.null(idx)) {
		td <- td[idx]
	}
	list(td = td, idx = idx)
}
