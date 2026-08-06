Partial output scheduling
=========================

``ProduceInputs`` declares two outputs: ``x.txt`` and ``y.txt``. ``UseX``
only consumes ``x.txt``, while ``UseY`` only consumes ``y.txt``.

To demonstrate scheduling only the output required by a named target, create
``x.txt`` without creating ``y.txt`` and turn off the conservative default::

    printf 'data for x\n' > x.txt
    gwf config set require_all_outputs false
    gwf run --dry-run UseX

The dry run reports ``Would submit UseX`` but does not report ``ProduceInputs``:
``x.txt`` is present and current, and ``y.txt`` is not needed for this request.
Running ``gwf run UseY`` instead, or ``gwf run UseX UseY``, schedules
``ProduceInputs`` because ``y.txt`` is required.

With the default setting, or after running::

    gwf config set require_all_outputs true

``gwf run --dry-run UseX`` also schedules ``ProduceInputs`` because its
declared ``y.txt`` output is missing.
