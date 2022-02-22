import ee


def accumulative_sum(image):
    bands = image.bandNames()

    def _sum(index, prev):
        curr_image = image.select([index]).unmask(0)
        prev_image = prev.get(-1).unmask(0)
        new_image = prev_image.add(curr_image)
        return prev.add(new_image)

    sums = bands.iterate(_sum, ee.List([ee.Image(0)]))
    return ee.ImageCollection.fromImages(sums).toBands().slice(1)


def find_stable_sequence(events, value, min_stability, max_masked, max_errors):
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

        test = (accumulative_matches
                .multiply(accumulative_errors.lte(max_errors))
                .reduce(ee.Reducer.max())
                .gte(min_stability))

        first_is_match = compare.select(0).eq(value)
        result = test.And(first_is_match)

        return result.selfMask()

    checks = ee.ImageCollection.fromImages(indices.map(_check)).toBands()

    matches = checks.multiply(indices_im)
    return matches.reduce(ee.Reducer.min())


def find_stable_sequences(events, values, min_stabilities, max_masked,
                          max_errors):
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
