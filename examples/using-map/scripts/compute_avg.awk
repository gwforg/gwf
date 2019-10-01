{
    sum += $1
    count += 1
} END {
    print sum / count
}