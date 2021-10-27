import java.text.SimpleDateFormat

Path path(filename) {
    file(filename, checkIfExists: true)
}

ArrayList<Map> read_tsv(Path path, List<String> names ) {
    path.toFile().readLines().with { lines ->
        lines.each {assert it.split('\t').size() == names.size() }
        lines.collect {
            [names, it.split('\t')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

ArrayList<Map> read_csv(Path path, List<String> names ) {
    path.toFile().readLines().with { lines ->
        lines.each {assert it.split(',').size() == names.size() }
        lines.collect {
            [names, it.split(',')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}
/* read_csv but with flexible column names and ordering */
ArrayList<Map> read_csv2(Path path, List<String> req_names, List<String> opt_names ) {
    path.toFile().readLines().with { lines ->
        names = lines[0].split(',') as ArrayList<String>
        assert req_names.all { names.contains(it) }
        // TODO - check opt_names
        lines.each {assert it.split(',').size() == names.size() }
        lines.collect {
            [names, it.split(',')].transpose().collectEntries { k, v -> [(k): v] }
        }
    }
}

ArrayList<ArrayList> get_families(ArrayList<Map> pedigree) {
    pedigree.groupBy { it.fid }
        .collect { k, v -> [
            k,
            v.findAll {it.phe == '2'}.collect {it.iid},
            v.findAll {it.phe == '1'}.collect {it.iid}
        ] }
}

String date_ymd() {
    date = new Date()
    sdf = new SimpleDateFormat("yyyy-MM-dd")
    sdf.format(date)
}