"""Methods that find the earliest occurence of a stable sequence in an image.

A stable sequence of a given value is at least X consecutive observations of
the value with at most N masked pixels and at most M errors occuring before the
Xth observation of the value of interest.
"""

import ee


def accumulative_sum(image):
    """Performs an accumulative sum on image treating masked pixels as 0s.

    Args:
        image: An ee.Image

    Returns:
        An ee.Image with the same number of bands as the input image where
        pixels values in band N represent the sum of all values in that pixel
        in the first N bands.
    """
    image = ee.Image(image)
    bands = image.bandNames()

    def _sum(index, prev):
        prev = ee.List(prev)
        curr_image = image.select([index]).unmask(0)
        prev_image = ee.Image(prev.get(-1)).unmask(0)
        new_image = prev_image.add(curr_image)
        return prev.add(new_image)

    sums = bands.iterate(_sum, ee.List([ee.Image(0)]))

    # slice(1) to drop the blank image iterate starts with
    return ee.ImageCollection.fromImages(sums).toBands().slice(1)


def find_stable_sequence(events, value, min_stability, max_masked, max_errors):
    """Finds the first stable sequence of value in events.

    Args:
        events: An ee.Image whose bands represent classifications through time.
        value: An integer to search for a stable sequence of.
        min_stability: An integer representing the minimum number of
            consecutive observations of value that must occur for a sequence
            to be considered stable.
        max_masked: An integer representing the maximum number of masked
            pixels that can occur in a sequence before min_stability
            observations and allow the sequence to still be considered stable.
        max_errors: An integer representing the maximum number of errors (i.e.
            pixels that do not equal value) that can occur in a sequence before
            min_stability observations and allwo the sequence to still be
            considered stable.

    Returns:
        An ee.Image with band whose values represent the band index in the
        input image where the earliest stable sequence of values occurs.
    """
    events = ee.Image(events)

    seq_length = ee.Number(max_masked).add(max_errors).add(min_stability)
    num_events = events.bandNames().length()
    indices = ee.List.sequence(0, num_events.subtract(1))
    indices_im = ee.Image.constant(indices)

    def _check(index):
        compare = events.slice(index, seq_length.add(index))
        matches = compare.eq(value)
        errors = compare.neq(value)

        accumulative_matches = accumulative_sum(matches)
        accumulative_errors = accumulative_sum(errors)

        # in order to search through masked pixels we search beyond
        # min_stability + max_errors elements, this can cause problems if there
        # are min_stability matches followed by max_errors + max_masked errors
        # which should match as it has a sequence of min_stability matches, but
        # will not because it has more than max_errors erros. To overcome this
        # we only consider how many matches occur before max_errors errors
        # occur.
        test = (accumulative_matches.multiply(
            accumulative_errors.lte(max_errors)).reduce(
                ee.Reducer.max()).gte(min_stability))

        # ensure that the first element matches to avoid off-by-one errors
        first_is_match = compare.select(0).eq(value)
        result = test.And(first_is_match)

        return result.selfMask()

    checks = ee.ImageCollection.fromImages(indices.map(_check)).toBands()

    matches = checks.multiply(indices_im)
    return matches.reduce(ee.Reducer.min())


def find_stable_sequences(events, values, min_stabilities, max_masked,
                          max_errors):
    """Find the first stable sequence of any of values in events.

    Args:
        events: An ee.Image whose bands represent classifications through time.
        values: An ee.List of integers representing each of the class values of
            interest to search for.
        min_stabilities: An ee.List of integers representing the stability for
            each of values.
        max_masked: An ee.List of integers representing the maximum allowable
            masked pixels in a stable sequence each of values.
        max_errors: An ee.List of intergers representing the maximum allowable
            errors in a stable sequence of each of values.

    Returns:
        An ee.Image with one band whose values represent the band index of the
        input image where the first stable sequence of values occurs.
    """
    events = ee.Image(events)
    values = ee.List(values)
    min_stabilities = ee.List(values)
    max_masked = ee.List(max_masked)
    max_errors = ee.List(max_errors)

    indices = ee.List.sequence(0, values.length().subtract(1))

    def _match(index):
        curr_value = values.get(index)
        curr_min_stability = min_stabilities.get(index)
        curr_max_masked = max_masked.get(index)
        curr_max_errors = max_errors.get(index)
        result = find_stable_sequence(events, curr_value, curr_min_stability,
                                      curr_max_masked, curr_max_errors)
        return result.int()

    matches = indices.map(_match).toBands()

    return matches.reduce(ee.Reducer.min(0))
