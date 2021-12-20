params.input = 'hello.txt'

condition = { it.readLines().size()>3 }
feedback_ch = Channel.create()
input_ch = Channel.fromPath(params.input).mix( feedback_ch.until(condition) )

process foo {
    input:
    file x from input_ch
    output:
    file 'foo.txt' into foo_ch
    script:
    """
    cat $x > foo.txt
    """
}

process bar {
    input:
    file x from foo_ch
    output:
    file 'bar.txt' into feedback_ch
    file 'bar.txt' into result_ch
    script:
    """
    cat $x > bar.txt
    echo World >> bar.txt
    """
}

result_ch.last().println { "Result:\n${it.text.indent(' ')}" }