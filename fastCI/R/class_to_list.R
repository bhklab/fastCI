class_to_list <-
function (x) 
{
    slot_names <- slotNames(x)
    l <- list()
    for (name in slot_names) {
        l[[name]] <- slot(x, name)
    }
    return(l)
}
