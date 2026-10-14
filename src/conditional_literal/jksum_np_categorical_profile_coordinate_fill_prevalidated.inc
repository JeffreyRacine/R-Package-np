static void np_categorical_profile_coordinate_fill_prevalidated(
  const NPCategoricalProfileKernelSpec *spec,
  int eval_index,
  int coordinate,
  int multiply,
  double *output)
{
  if(coordinate < spec->nunordered) {
    const int kernel = spec->kernel_unordered[coordinate] +
      OP_UFUN_OFFSETS[spec->operator_code[coordinate]];

    np_ukernelv(kernel,
                spec->train_unordered[coordinate],
                spec->ntrain,
                multiply,
                spec->eval_unordered[coordinate][eval_index],
                spec->lambda[coordinate],
                spec->num_categories[coordinate],
                output,
                NULL,
                NP_CF316_SKIP_UPPER(0));
  } else {
    const int ordered_coordinate = coordinate - spec->nunordered;
    const int base_kernel = spec->kernel_ordered[ordered_coordinate];
    const int operator_code = spec->operator_code[coordinate];

    /* Preserve the incumbent ordered Li--Racine convolution owner. */
    if(operator_code == OP_CONVOLUTION && base_kernel == 1) {
      np_convol_okernelv(base_kernel,
                         spec->train_ordered[ordered_coordinate],
                         spec->ntrain,
                         multiply,
                         spec->eval_ordered[ordered_coordinate][eval_index],
                         spec->lambda[coordinate],
                         spec->num_categories[coordinate],
                         spec->category_values[coordinate],
                         output,
                         0);
    } else {
      const int kernel = base_kernel + OP_OFUN_OFFSETS[operator_code];

      np_okernelv(kernel,
                  spec->train_ordered[ordered_coordinate],
                  spec->ntrain,
                  multiply,
                  spec->eval_ordered[ordered_coordinate][eval_index],
                  spec->lambda[coordinate],
                  spec->category_values[coordinate],
                  spec->num_categories[coordinate],
                  output,
                  NULL,
                  0);
    }
  }
}
