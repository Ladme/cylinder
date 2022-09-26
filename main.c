// Released under MIT License.
// Copyright (c) 2022 Ladislav Bartos

#include <unistd.h>
#include <groan.h>

// frequency of printing during the calculation
static const int PROGRESS_FREQ = 10000;
// version of the code
static const char VERSION[] = "v2022/09/26";

static const float PI = 3.141592f;

/*
 * Parses command line arguments.
 * Returns zero, if parsing has been successful. Else returns non-zero.
 */
int get_arguments(
        int argc, 
        char **argv,
        char **gro_file,
        char **xtc_file,
        char **ndx_file,
        char **output_file,
        char **selection,
        char **position,
        dimension_t *axis,
        float *radius,
        float *height,
        int *dx) 
{
    int gro_specified = 0, selection_specified = 0, position_specified = 0, axis_specified = 0;

    int opt = 0;
    while((opt = getopt(argc, argv, "c:f:n:o:s:p:xyzr:e:d:h")) != -1) {
        switch (opt) {
        // help
        case 'h':
            return 1;
        // gro file to read
        case 'c':
            *gro_file = optarg;
            gro_specified = 1;
            break;
        // xtc file to read
        case 'f':
            *xtc_file = optarg;
            break;
        // ndx file to read
        case 'n':
            *ndx_file = optarg;
            break;
        // output file
        case 'o':
            *output_file = optarg;
            break;
        // atom selection
        case 's':
            *selection = optarg;
            selection_specified = 1;
            break;
        // selection(s) defining position (center) of the cylinder
        case 'p':
            *position = optarg;
            position_specified = 1;
            break;
        case 'x':
        case 'y':
        case 'z':
            if (axis_specified) {
                fprintf(stderr, "Multiple axes specified.\n");
                return 1;
            }

            axis_specified = 1;
            if (opt == 'x') *axis = x;
            else if (opt == 'y') *axis = y;
            else *axis = z;
            break;
        // radius for the cylinder
        case 'r':
            *radius = atof(optarg);
            if (*radius <= 0) {
                fprintf(stderr, "Cylinder radius must be >0, not %f.\n", *radius);
                return 1;
            }
            break;
        // height of the cylinder
        case 'e':
            *height = atof(optarg);
            if (*height <= 0) {
                fprintf(stderr, "Cylinder height must be >0, not %f.\n", *height);
                return 1;
            }
            break;
        // spacing along the axis
        case 'd':
            *dx = atoi(optarg);
            if (*dx <= 0) {
                fprintf(stderr, "Spacing must be > 0, not %d.\n", *dx);
                return 1;
            }
            break;
        default:
            //fprintf(stderr, "Unknown command line option: %c.\n", opt);
            return 1;
        }
    }

    if (!gro_specified || !selection_specified || !position_specified) {
        fprintf(stderr, "Gro file, selection and position specification must always be supplied.\n");
        return 1;
    }
    return 0;
}

void print_usage(const char *program_name)
{
    printf("Usage: %s -c GRO_FILE -s SELECTION -p CYLINDER_POSITION [OPTION]...\n", program_name);
    printf("\nOPTIONS\n");
    printf("-h               print this message and exit\n");
    printf("-c STRING        gro file to read\n");
    printf("-f STRING        xtc file to read (optional)\n");
    printf("-n STRING        ndx file to read (optional, default: index.ndx)\n");
    printf("-o STRING        output file (default: cylinder.xvg)\n");
    printf("-s STRING        selection of atoms to analyze\n");
    printf("-p STRING        selection of atoms defining cylinder position\n");
    printf("-x/y/z           direction of the main axis of the cylinder (default: z)\n");
    printf("-r FLOAT         radius of the cylinder in nm (default: 2.5)\n");
    printf("-e FLOAT         height of the cylinder in nm (default: height of the box in gro file)\n");
    printf("-d INTEGER       grid spacing along the specified axis in points per nm (default: 10)\n");
    printf("\n");
}

/*
 * Prints parameters that the program will use for the cylinder density calculation.
 */
void print_arguments(
        const char *gro_file, 
        const char *xtc_file,
        const char *ndx_file,
        const char *output_file,
        const char *selection,
        char **parsed,
        const int n_parsed,
        const dimension_t axis,
        const float radius,
        const float height,
        const int dx)
{
    printf("Parameters for Cylinder calculation:\n");
    printf(">>> gro file:         %s\n", gro_file);
    if (xtc_file != NULL) printf(">>> xtc file:         %s\n", xtc_file);
    printf(">>> ndx file:         %s\n", ndx_file);
    printf(">>> output file:      %s\n", output_file);
    printf(">>> selection:        %s\n", selection);
    
    if (n_parsed == 1) {
        printf(">>> position:         %s\n", parsed[0]);
    }
    else {
        strstrip(parsed[0]);
        strstrip(parsed[1]);
        strstrip(parsed[2]);
        printf(">>> x-position:       %s\n", parsed[0]);
        printf(">>> y-position:       %s\n", parsed[1]);
        printf(">>> z-position:       %s\n", parsed[2]);
    }

    if (axis == x) printf(">>> cylinder axis:    x\n");
    else if (axis == y) printf(">>> cylinder axis:    y\n");
    else printf(">>> cylinder axis:    z\n");

    printf(">>> cylinder radius:  %f\n", radius);
    printf(">>> cylinder height:  %f\n", height);
    printf(">>> spacing:          %d\n", dx);
    printf("\n");
}

/*
 * Converts index of an array to coordinate.
 */
static inline float index2coor(int x, float minx, int dx)
{
    return (float) x / dx + minx;
}
        
/*
 * Converts coordinate to an index array.
 */
static inline size_t coor2index(float x, float minx, int dx)
{
    return (size_t) roundf((x - minx) * dx);
}


/*! @brief Calculates position of the cylinder. */
void calc_position(
        system_t *system, 
        vec_t cyl_pos,
        const atom_selection_t *positionx, 
        const atom_selection_t *positiony,
        const atom_selection_t *positionz,
        const int equal_xy,
        const int equal_xz,
        const int equal_yz,
        const int equal_xyz)
{   
    if (equal_xyz) {
        center_of_geometry(positionx, cyl_pos, system->box);

    } else if (equal_xy) {
        vec_t center_z = {0.0};
        center_of_geometry(positionx, cyl_pos, system->box);
        center_of_geometry(positionz, center_z, system->box);
        
        cyl_pos[2] = center_z[2];

    } else if (equal_xz) {
        vec_t center_y = {0.0};
        center_of_geometry(positionx, cyl_pos, system->box);
        center_of_geometry(positiony, center_y, system->box);

        cyl_pos[1] = center_y[1];

    } else if (equal_yz) {
        vec_t center_x = {0.0};
        center_of_geometry(positiony, cyl_pos, system->box);
        center_of_geometry(positionx, center_x, system->box);

        cyl_pos[0] = center_x[0];

    } else {
        vec_t center_x = {0.0};
        vec_t center_y = {0.0};
        vec_t center_z = {0.0};

        center_of_geometry(positionx, center_x, system->box);
        center_of_geometry(positiony, center_y, system->box);
        center_of_geometry(positionz, center_z, system->box);

        cyl_pos[0] = center_x[0];
        cyl_pos[1] = center_y[1];
        cyl_pos[2] = center_z[2];
    }
}

/*! Calculates the number of selected atoms inside the specified cylinder along the specified axis. */
void calc_density(
        system_t *system, 
        size_t *density,
        const atom_selection_t *selection,
        const vec_t cylinder_position,
        const plane_t cylinder_plane,
        const dimension_t cylinder_axis,
        const float radius,
        const float half_height,
        const int dx,
        const size_t array_size)
{  
    // assign atoms to slices
    for (size_t i = 0; i < selection->n_atoms; ++i) {
        atom_t *atom = selection->atoms[i];

        // filter out atoms that are not inside of the cylinder
        if (distance2D(atom->position, cylinder_position, cylinder_plane, system->box) > radius) continue;
        float distance = distance1D(atom->position, cylinder_position, cylinder_axis, system->box);
        if (fabsf(distance) > half_height) continue;

        size_t index = coor2index(distance, (0.5 / dx) - half_height, dx);
        if (index >= array_size) continue;

        //printf("Distance: %f, Index: %lu, Distance index: %f\n", distance, index, index2coor(index, (0.5 / dx) - half_height, dx));
        // increase the number of atoms in the corresponding slice
        density[index]++;
    }
}

int main(int argc, char **argv)
{
    int return_code = 0;

    printf("\n");
    // get command line arguments
    char *gro_file = NULL;
    char *xtc_file = NULL;
    char *ndx_file = "index.ndx";
    char *output_file = "cylinder.xvg";
    char *selection = NULL;
    char *position = NULL;
    dimension_t axis = z;
    float radius = 2.5;
    float height = 0.0;
    int dx = 10;

    if (get_arguments(argc, argv, &gro_file, &xtc_file, &ndx_file, &output_file, &selection, &position, &axis, &radius, &height, &dx) != 0) {
        print_usage(argv[0]);
        return 1;
    }

    // read gro file
    system_t *system = load_gro(gro_file);
    if (system == NULL) return 1;

    // select all atoms
    atom_selection_t *all = select_system(system);

    // read ndx file
    dict_t *ndx_groups = read_ndx(ndx_file, system);

    // select atoms
    atom_selection_t *atoms = smart_select(all, selection, ndx_groups);
    if (atoms == NULL || atoms->n_atoms == 0) {
        fprintf(stderr, "No atoms ('%s') found.\n", selection);
        dict_destroy(ndx_groups);
        free(all);
        free(atoms);
        free(system);
        return 1;
    }

    // parse position specifiers
    char **parsed = NULL;
    char *position_original = strdup(position);
    size_t n_items = strsplit(position, &parsed, ";");
    atom_selection_t *positionx = NULL;
    atom_selection_t *positiony = NULL;
    atom_selection_t *positionz = NULL;
    int equal_xy = 0;
    int equal_xz = 0;
    int equal_yz = 0;
    int equal_xyz = 0;

    // one specification for all dimensions
    if (n_items == 1) {
        positionx = smart_select(all, position_original, ndx_groups);
        if (positionx == NULL || positionx->n_atoms == 0) {
            fprintf(stderr, "No atoms ('%s') found.\n", position_original);
            goto position_fail;
        }

        equal_xyz = 1;
    // specifications for individual dimensions
    } else if (n_items == 3) {
        positionx = smart_select(all, parsed[0], ndx_groups);
        if (positionx == NULL || positionx->n_atoms == 0) {
            fprintf(stderr, "No atoms ('%s') found.\n", parsed[0]);
            goto position_fail;
        }

        positiony = smart_select(all, parsed[1], ndx_groups);
        if (positiony == NULL || positiony->n_atoms == 0) {
            fprintf(stderr, "No atoms ('%s') found.\n", parsed[1]);
            goto position_fail;
        }

        positionz = smart_select(all, parsed[2], ndx_groups);
        if (positionz == NULL || positionz->n_atoms == 0) {
            fprintf(stderr, "No atoms ('%s') found.\n", parsed[2]);
            goto position_fail;
        }

        if (selection_compare(positionx, positiony)) equal_xy = 1;
        if (selection_compare(positionx, positionz)) equal_xz = 1;
        if (selection_compare(positiony, positionz)) equal_yz = 1;

        if (equal_xy && equal_xz) {
            equal_xyz = 1;
            strstrip(parsed[0]);
            printf("Note: You can use '%s' instead of '%s' in the position specification.\n\n", parsed[0], position_original);
        }

    // different number results in an error
    } else {
        fprintf(stderr, "Could not parse selection '%s'\n", position_original);

        position_fail:
        dict_destroy(ndx_groups);
        free(all);
        free(atoms);
        free(system);
        free(positionx);
        free(positiony);
        free(positionz);
        free(parsed);
        free(position_original);
        return 1;
    }

    free(position_original);

    // get default height of the cylinder, if not supplied
    if (height == 0) {
        if (axis == x) height = system->box[0];
        else if (axis == y) height = system->box[1];
        else height = system->box[2];
    }
    float half_height = height / 2;

    print_arguments(gro_file, xtc_file, ndx_file, output_file, selection, parsed, n_items, axis, radius, height, dx);

    free(parsed);

    // prepare density array
    size_t array_size = (size_t) roundf( height * dx );
    size_t *density = calloc(array_size, sizeof(size_t));

    // get plane of the cylinder
    plane_t plane = xy;
    if (axis == x) plane = yz;
    else if (axis == y) plane = xz;

    // open output file & prepare header
    FILE *output = fopen(output_file, "w");
    if (output == NULL) {
        fprintf(stderr, "File %s could not be opened.\n", output_file);
        return_code = 1;
        goto program_end;
    }
    
    fprintf(output, "# Generated with cylinder (C Cylinder Density Calculator) %s\n", VERSION);
    fprintf(output, "# Command line: ");
    for (int i = 0; i < argc; ++i) {
        fprintf(output, "%s ", argv[i]);
    }
    fprintf(output, "\n");
    fprintf(output, "# See the average number of atoms inside the entire cylinder at the end of this file.\n");

    fprintf(output, "@    title \"Density of atoms in a cylinder\"\n");
    if (axis == x) fprintf(output, "@    xaxis label \"x-axis\"\n");
    else if (axis == y) fprintf(output, "@    xaxis label \"y-axis\"\n");
    else fprintf(output, "@    xaxis label \"z-axis\"\n");
    fprintf(output, "@    yaxis label \"density [particles per nm^3]\"\n");
    fprintf(output, "@    s1 legend \"%s\"\n", selection);

    size_t n_frames = 0;
    // in case the xtc file is not supplied, analyze only the gro file
    if (xtc_file == NULL) {
        vec_t cyl_pos = {0.0};
        calc_position(system, cyl_pos, positionx, positiony, positionz, equal_xy, equal_xz, equal_yz, equal_xyz);
        calc_density(system, density, atoms, cyl_pos, plane, axis, radius, half_height, dx, array_size);
        ++n_frames;
    // open and read xtc file
    } else {
        XDRFILE *xtc = xdrfile_open(xtc_file, "r");

        if (xtc == NULL) {
            fprintf(stderr, "File %s could not be read as an xtc file.\n", xtc_file);
            return_code = 1;
            goto program_end;
        }

        if (!validate_xtc(xtc_file, (int) system->n_atoms)) {
            fprintf(stderr, "Number of atoms in %s does not match %s.\n", xtc_file, gro_file);
            return_code = 1;
            xdrfile_close(xtc);
            goto program_end;
        }
        
        while (read_xtc_step(xtc, system) == 0) {
            
            // print info about the progress of reading
            if ((int) system->time % PROGRESS_FREQ == 0) {
                printf("Step: %d. Time: %.0f ps\r", system->step, system->time);
                fflush(stdout);
            }
            
            vec_t cyl_pos = {0.0};
            calc_position(system, cyl_pos, positionx, positiony, positionz, equal_xy, equal_xz, equal_yz, equal_xyz);
            calc_density(system, density, atoms, cyl_pos, plane, axis, radius, half_height, dx, array_size);
            ++n_frames;
        }

        xdrfile_close(xtc);
    }

    // write output
    size_t total_number = 0;
    float slice_volume = (1.0 / dx) * PI * radius * radius;
    for (size_t i = 0; i < array_size; ++i) {
        total_number += density[i];
        fprintf(output, "%f %f\n", index2coor(i, (0.5 / dx) - half_height, dx), (float) density[i] / n_frames / slice_volume);
    }

    fprintf(output, "# Average number of particles inside the cylinder: %f\n", (float) total_number / n_frames);

    printf("Output written into file %s.\n", output_file);

    program_end:
    if (output != NULL) fclose(output);
    free(density);
    dict_destroy(ndx_groups);
    free(all);
    free(atoms);
    free(positionx);
    free(positiony);
    free(positionz);
    free(system);

    return return_code;
}