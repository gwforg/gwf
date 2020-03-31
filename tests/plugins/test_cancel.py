from gwf.cli import main


def test_cancel_one_target(cli_runner):
    result = cli_runner.invoke(
        main,
        [
            "-b",
            "testing",
            "-f",
            "examples/simple-workflow/workflow.py",
            "cancel",
            "Target1",
        ],
    )
    assert result.output == "Cancelling target Target1\n"


def test_cancel_two_targets(cli_runner):
    result = cli_runner.invoke(
        main,
        [
            "-b",
            "testing",
            "-f",
            "examples/simple-workflow/workflow.py",
            "cancel",
            "Target1",
            "Target2",
        ],
    )
    lines = result.output.split("\n")
    assert len(lines) == 3
    assert "Cancelling target Target1" in lines
    assert "Cancelling target Target2" in lines


def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(
    cli_runner,
):
    result = cli_runner.invoke(
        main,
        ["-b", "testing", "-f", "examples/simple-workflow/workflow.py", "cancel"],
        input="y",
    )
    lines = result.output.split("\n")
    assert len(lines) == 4
    assert "Cancelling target Target1" in lines
    assert "Cancelling target Target2" in lines


def test_cancel_no_targets_specified_should_ask_for_confirmation_and_abort_if_not_approved(
    cli_runner,
):
    result = cli_runner.invoke(
        main,
        ["-b", "testing", "-f", "examples/simple-workflow/workflow.py", "cancel"],
        input="N",
    )
    assert "Aborted!\n" in result.output
