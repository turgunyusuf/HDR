/************************************************************************/
/*                                                                      */
/*               Copyright 2001-2002 by Gunnar Kedenburg                */
/*        Copyright 2012 Christoph Spiel and Ullrich Koethe             */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/


/*!
 * \file  impex.hxx
 * \brief image import and export functions
 *
 * This module implements functions importImage() and exportImage().
 * The matching implementation for any given datatype is selected by
 * template meta code.
 *
 */

#ifndef VIGRA_IMPEX_HXX
#define VIGRA_IMPEX_HXX

#include "stdimage.hxx"
#include "imageinfo.hxx"
#include "impexbase.hxx"

namespace vigra
{
/** \addtogroup VigraImpex
 * @{
*/
    namespace detail
    {
        template <class ValueType,
                  class ImageIterator, class ImageAccessor>
        void
        read_image_band(Decoder* decoder,
                        ImageIterator image_iterator, ImageAccessor image_accessor)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;

            const unsigned width(decoder->getWidth());
            const unsigned height(decoder->getHeight());
            const unsigned offset(decoder->getOffset());

            for (unsigned y = 0U; y != height; ++y)
            {
                decoder->nextScanline();

                const ValueType* scanline = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));

                ImageRowIterator is(image_iterator.rowIterator());
                const ImageRowIterator is_end(is + width);

                while (is != is_end)
                {
                    image_accessor.set(*scanline, is);
                    scanline += offset;
                    ++is;
                }

                ++image_iterator.y;
            }
        }


        template <class ValueType,
                  class ImageIterator, class ImageAccessor>
        void
        read_image_bands(Decoder* decoder,
                         ImageIterator image_iterator, ImageAccessor image_accessor)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;

            const unsigned width(decoder->getWidth());
            const unsigned height(decoder->getHeight());
            const unsigned offset(decoder->getOffset());
            const unsigned accessor_size(image_accessor.size(image_iterator));

            // OPTIMIZATION: Specialization for the most common case
            // of an RGB-image, i.e. 3 channels.
            if (accessor_size == 3U)
            {
                const ValueType* scanline_0;
                const ValueType* scanline_1;
                const ValueType* scanline_2;

                for (unsigned y = 0U; y != height; ++y)
                {
                    decoder->nextScanline();

                    scanline_0 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<const ValueType*>(decoder->currentScanlineOfBand(2));

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        image_accessor.setComponent(*scanline_0, is, 0);
                        image_accessor.setComponent(*scanline_1, is, 1);
                        image_accessor.setComponent(*scanline_2, is, 2);

                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;

                        ++is;
                    }

                    ++image_iterator.y;
                }
            }
            else
            {
                std::vector<const ValueType*> scanlines(accessor_size);

                for (unsigned y = 0U; y != height; ++y)
                {
                    decoder->nextScanline();

                    for (unsigned i = 0U; i != accessor_size; ++i)
                    {
                        scanlines[i] = static_cast<const ValueType*>(decoder->currentScanlineOfBand(i));
                    }

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        for (unsigned i = 0U; i != accessor_size; ++i)
                        {
                            image_accessor.setComponent(*scanlines[i], is, static_cast<int>(i));
                            scanlines[i] += offset;
                        }
                        ++is;
                    }

                    ++image_iterator.y;
                }
            }
        }


        template <class ImageIterator, class ImageAccessor>
        void
        importImage(const ImageImportInfo& import_info,
                    ImageIterator image_iterator, ImageAccessor image_accessor,
                    /* isScalar? */ VigraTrueType)
        {
            VIGRA_UNIQUE_PTR<Decoder> decoder(vigra::decoder(import_info));

            switch (pixel_t_of_string(decoder->getPixelType()))
            {
            case UNSIGNED_INT_8:
                read_image_band<UInt8>(decoder.get(), image_iterator, image_accessor);
                break;
            case UNSIGNED_INT_16:
                read_image_band<UInt16>(decoder.get(), image_iterator, image_accessor);
                break;
            case UNSIGNED_INT_32:
                read_image_band<UInt32>(decoder.get(), image_iterator, image_accessor);
                break;
            case SIGNED_INT_16:
                read_image_band<Int16>(decoder.get(), image_iterator, image_accessor);
                break;
            case SIGNED_INT_32:
                read_image_band<Int32>(decoder.get(), image_iterator, image_accessor);
                break;
            case IEEE_FLOAT_32:
                read_image_band<float>(decoder.get(), image_iterator, image_accessor);
                break;
            case IEEE_FLOAT_64:
                read_image_band<double>(decoder.get(), image_iterator, image_accessor);
                break;
            default:
                vigra_fail("detail::importImage<scalar>: not reached");
            }

            decoder->close();
        }


        template <class ImageIterator, class ImageAccessor>
        void
        importImage(const ImageImportInfo& import_info,
                    ImageIterator image_iterator, ImageAccessor image_accessor,
                    /* isScalar? */ VigraFalseType)
        {
            VIGRA_UNIQUE_PTR<Decoder> decoder(vigra::decoder(import_info));

            switch (pixel_t_of_string(decoder->getPixelType()))
            {
            case UNSIGNED_INT_8:
                read_image_bands<UInt8>(decoder.get(), image_iterator, image_accessor);
                break;
            case UNSIGNED_INT_16:
                read_image_bands<UInt16>(decoder.get(), image_iterator, image_accessor);
                break;
            case UNSIGNED_INT_32:
                read_image_bands<UInt32>(decoder.get(), image_iterator, image_accessor);
                break;
            case SIGNED_INT_16:
                read_image_bands<Int16>(decoder.get(), image_iterator, image_accessor);
                break;
            case SIGNED_INT_32:
                read_image_bands<Int32>(decoder.get(), image_iterator, image_accessor);
                break;
            case IEEE_FLOAT_32:
                read_image_bands<float>(decoder.get(), image_iterator, image_accessor);
                break;
            case IEEE_FLOAT_64:
                read_image_bands<double>(decoder.get(), image_iterator, image_accessor);
                break;
            default:
                vigra_fail("vigra::detail::importImage<non-scalar>: not reached");
            }

            decoder->close();
        }

        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler>
        void
        write_image_band(Encoder* encoder,
                         ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                         const ImageScaler& image_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef typename ImageAccessor::value_type ImageValueType;

            typedef RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_band: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_band: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(1);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default
            // constructor to allow classes ImageIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);

            for (unsigned y = 0U; y != height; ++y)
            {
                ValueType* scanline = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));

                ImageRowIterator is(image_iterator.rowIterator());
                const ImageRowIterator is_end(is + width);

                while (is != is_end)
                {
                    *scanline = explicit_cast::cast(image_scaler(image_accessor(is)));
                    scanline += offset;
                    ++is;
                }

                encoder->nextScanline();

                ++image_iterator.y;
            }
        }


        template<class ValueType,
                 class ImageIterator, class ImageAccessor, class ImageScaler>
        void
        write_image_bands(Encoder* encoder,
                          ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                          const ImageScaler& image_scaler)
        {
            typedef typename ImageIterator::row_iterator ImageRowIterator;
            typedef RequiresExplicitCast<ValueType> explicit_cast;

            vigra_precondition(image_lower_right.x >= image_upper_left.x,
                               "vigra::detail::write_image_bands: negative width");
            vigra_precondition(image_lower_right.y >= image_upper_left.y,
                               "vigra::detail::write_image_bands: negative height");

            const unsigned width(static_cast<unsigned>(image_lower_right.x - image_upper_left.x));
            const unsigned height(static_cast<unsigned>(image_lower_right.y - image_upper_left.y));
            const unsigned accessor_size(image_accessor.size(image_upper_left));

            encoder->setWidth(width);
            encoder->setHeight(height);
            encoder->setNumBands(accessor_size);
            encoder->finalizeSettings();

            const unsigned offset(encoder->getOffset()); // correct offset only _after_ finalizeSettings()

            // IMPLEMENTATION NOTE: We avoid calling the default
            // constructor to allow classes ImageIterator that do not
            // define one.
            ImageIterator image_iterator(image_upper_left);

            // OPTIMIZATION: Specialization for the most common case
            // of an RGB-image, i.e. 3 channels.
            if (accessor_size == 3U)
            {
                ValueType* scanline_0;
                ValueType* scanline_1;
                ValueType* scanline_2;

                for (unsigned y = 0U; y != height; ++y)
                {
                    scanline_0 = static_cast<ValueType*>(encoder->currentScanlineOfBand(0));
                    scanline_1 = static_cast<ValueType*>(encoder->currentScanlineOfBand(1));
                    scanline_2 = static_cast<ValueType*>(encoder->currentScanlineOfBand(2));

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        *scanline_0 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 0)));
                        *scanline_1 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 1)));
                        *scanline_2 = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, 2)));

                        scanline_0 += offset;
                        scanline_1 += offset;
                        scanline_2 += offset;

                        ++is;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                }
            }
            else
            {
                std::vector<ValueType*> scanlines(accessor_size);

                for (unsigned y = 0U; y != height; ++y)
                {
                    for (unsigned i = 0U; i != accessor_size; ++i)
                    {
                        scanlines[i] = static_cast<ValueType*>(encoder->currentScanlineOfBand(i));
                    }

                    ImageRowIterator is(image_iterator.rowIterator());
                    const ImageRowIterator is_end(is + width);

                    while (is != is_end)
                    {
                        for (unsigned i = 0U; i != accessor_size; ++i)
                        {
                            *scanlines[i] = explicit_cast::cast(image_scaler(image_accessor.getComponent(is, static_cast<int>(i))));
                            scanlines[i] += offset;
                        }
                        ++is;
                    }

                    encoder->nextScanline();

                    ++image_iterator.y;
                }
            }
        }


        template <class ImageIterator, class ImageAccessor>
        void
        exportImage(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                    const ImageExportInfo& export_info,
                    /* isScalar? */ VigraTrueType)
        {
            typedef typename ImageAccessor::value_type ImageValueType;

            VIGRA_UNIQUE_PTR<Encoder> encoder(vigra::encoder(export_info));

            std::string pixel_type(export_info.getPixelType());
            const bool downcast(negotiatePixelType(encoder->getFileType(), TypeAsString<ImageValueType>::result(), pixel_type));
            const pixel_t type(pixel_t_of_string(pixel_type));

            encoder->setPixelType(pixel_type);

            const range_t image_source_range(find_source_value_range(export_info,
                                                                     image_upper_left, image_lower_right, image_accessor));
            const range_t destination_range(find_destination_value_range(export_info, type));

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_band<UInt8>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_band<UInt16>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_band<UInt32>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_band<Int16>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_band<Int32>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_band<float>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_band<double>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_band<UInt8>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_band<UInt16>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_band<UInt32>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_16:
                    write_image_band<Int16>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_32:
                    write_image_band<Int32>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_band<float>(encoder.get(),
                                            image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_band<double>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<scalar>: not reached");
                }
            }

            encoder->close();
        }


        template <class ImageIterator, class ImageAccessor>
        void
        exportImage(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                    const ImageExportInfo& export_info,
                    /* isScalar? */ VigraFalseType)
        {
            typedef typename ImageAccessor::value_type ImageBaseType;
            typedef typename ImageBaseType::value_type ImageValueType;

            VIGRA_UNIQUE_PTR<Encoder> encoder(vigra::encoder(export_info));

            std::string pixel_type(export_info.getPixelType());
            const bool downcast(negotiatePixelType(encoder->getFileType(), TypeAsString<ImageValueType>::result(), pixel_type));
            const pixel_t type(pixel_t_of_string(pixel_type));

            encoder->setPixelType(pixel_type);

            vigra_precondition(isBandNumberSupported(encoder->getFileType(), image_accessor.size(image_upper_left)),
                               "exportImage(): file format does not support requested number of bands (color channels)");

            const range_t image_source_range(find_source_value_range(export_info,
                                                                     image_upper_left, image_lower_right, image_accessor));
            const range_t destination_range(find_destination_value_range(export_info, pixel_t_of_string(pixel_type)));

            if ((downcast || export_info.hasForcedRangeMapping()) &&
                (image_source_range.first != destination_range.first || image_source_range.second != destination_range.second))
            {
                const linear_transform image_rescaler(image_source_range, destination_range);

                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<UInt8>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<UInt16>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<UInt32>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_16:
                    write_image_bands<Int16>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case SIGNED_INT_32:
                    write_image_bands<Int32>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, image_rescaler);
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<non-scalar>: not reached");
                }
            }
            else
            {
                switch (type)
                {
                case UNSIGNED_INT_8:
                    write_image_bands<UInt8>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_16:
                    write_image_bands<UInt16>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case UNSIGNED_INT_32:
                    write_image_bands<UInt32>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_16:
                    write_image_bands<Int16>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case SIGNED_INT_32:
                    write_image_bands<Int32>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_32:
                    write_image_bands<float>(encoder.get(),
                                             image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                case IEEE_FLOAT_64:
                    write_image_bands<double>(encoder.get(),
                                              image_upper_left, image_lower_right, image_accessor, identity());
                    break;
                default:
                    vigra_fail("vigra::detail::exportImage<non-scalar>: not reached");
                }
            }

            encoder->close();
        }
    }  // end namespace detail

    /*!
     * \brief Read the image specified by the given \ref
     * vigra::ImageImportInfo object.
     *
     * <B>Declarations</B>
     *
     * Pass arguments explicitly:
     * \code
     * namespace vigra {
     *     template <class ImageIterator, class Accessor>
     *     void
     *     importImage(const ImageImportInfo& importInfo,
     *                 ImageIterator imageIterator, Accessor imageAccessor)
     * }
     * \endcode
     *
     * Use argument objects in conjunction with \ref ArgumentObjectFactories :
     * \code
     * namespace vigra {
     *     template <class ImageIterator, class Accessor>
     *     inline void
     *     importImage(const ImageImportInfo& importInfo,
     *                 const pair<ImageIterator, Accessor>& image)
     * }
     * \endcode
     *
     * <B>Usage</B>
     *
     * <B>\#include \<vigra/impex.hxx\></B>
     *
     * Namespace: vigra
     *
     * \code
     *     ImageImportInfo info("myimage.gif");
     *
     *     if (info.isGrayscale())
     *     {
     *         // create byte image of appropriate size
     *         BImage image(info.width(), info.height());
     *
     *         importImage(info, destImage(image));
     *         ...
     *     }
     *     else
     *     {
     *         // create byte RGB image of appropriate size
     *         BRGBImage image(info.width(), info.height());
     *
     *         importImage(info, destImage(image));
     *         ...
     *     }
     * \endcode
     *
     * <B>Preconditions</B>
     *
     * - The image file must be readable and
     * - the file type must be one of the following:
     *
     * | Type | Extension | Name                                                       | Support Library |
     * |------|-----------|------------------------------------------------------------|-----------------|
     * | BMP  | bmp       | Microsoft Windows bitmap image file                        | |
     * | EXR  | exr       | OpenEXR high dynamic range image format                    | libopenexr |
     * | GIF  | gif       | CompuServe graphics interchange format, 8-bit color        | |
     * | HDR  | hdr       | Radiance RGBE high dynamic range image format              | libexr? |
     * | JPEG | jpg, jpeg | Joint Photographic Experts Group JFIF format, 24-bit color | libjpeg |
     * | PBM  | pbm       | Portable bitmap format (black and white)                   | |
     * | PGM  | pgm       | Portable graymap format (gray scale)                       | |
     * | PNG  | png       | Portable Network Graphic                                   | libpng |
     * | PNM  | pnm       | Portable anymap                                            | |
     * | PPM  | ppm       | Portable pixmap format (color)                             | |
     * | SUN  | ras       | SUN Rasterfile                                             | |
     * | TIFF | tif, tiff | Tagged Image File Format                                   | libtiff |
     * | VIFF | xv        | Khoros Visualization image file                            | |
     */
    doxygen_overloaded_function(template <...> inline void importImage)


    template <class ImageIterator, class ImageAccessor>
    inline void
    importImage(const ImageImportInfo& import_info,
                ImageIterator image_iterator, ImageAccessor image_accessor)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename NumericTraits<ImageValueType>::isScalar is_scalar;

        detail::importImage(import_info,
                    image_iterator, image_accessor,
                    is_scalar());
    }


    template <class ImageIterator, class ImageAccessor>
    inline void
    importImage(const ImageImportInfo& import_info,
                const vigra::pair<ImageIterator, ImageAccessor>& image)
    {
        importImage(import_info,
                    image.first, image.second);
    }

    /*!
     * \brief Write an image given a \ref vigra::ImageExportInfo object.
     *
     * If the file format to be exported to supports the pixel type of
     * the source image, the pixel type will be kept
     * (e.g. <tt>float</tt> can be stored as TIFF without conversion,
     * in contrast to most other image export toolkits).  Otherwise,
     * the pixel values are transformed to the range 0..255 and
     * converted to <tt>unsigned char</tt>.  Currently, the following
     * file formats are supported.  The pixel types given in brackets
     * are those that are written without conversion:
     *     - BMP: Microsoft Windows bitmap image file (pixel types: UINT8 as gray and RGB);
     *     - GIF: CompuServe graphics interchange format, 8-bit color (pixel types: UINT8 as gray and RGB);
     *     - JPEG: Joint Photographic Experts Group JFIF format, compressed 24-bit color
     *             (pixel types: UINT8 as gray and RGB), only available if libjpeg is installed;
     *     - PNG: Portable Network Graphic (pixel types: UINT8 and UINT16 with up to 4 channels),
     *             only available if libpng is installed;
     *     - PBM: Portable bitmap format (black and white);
     *     - PGM: Portable graymap format (pixel types: UINT8, INT16, INT32 as gray scale);
     *     - PNM: Portable anymap (pixel types: UINT8, INT16, INT32 as gray and RGB);
     *     - PPM: Portable pixmap format (pixel types: UINT8, INT16, INT32 as RGB);
     *     - SUN: SUN Rasterfile (pixel types: UINT8 as gray and RGB);
     *     - TIFF: Tagged Image File Format
     *           (pixel types: UINT8, INT16, INT32, FLOAT, DOUBLE with up to 4 channels),
     *           only available if libtiff is installed;
     *     - VIFF: Khoros Visualization image file
     *           (pixel types: UINT8, INT16, INT32, FLOAT, DOUBLE with arbitrary many channels);
     *
     * <B>Declarations</B>
     *
     * Pass arguments explicitly:
     * \code
     * namespace vigra {
     *     template <class ImageIterator, class ImageAccessor>
     *     void
     *     exportImage(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
     *                 const ImageExportInfo& exportInfo)
     * }
     * \endcode
     *
     * Use argument objects in conjunction with \ref ArgumentObjectFactories :
     * \code
     *     namespace vigra {
     *         template <class ImageIterator, class ImageAccessor>
     *         void exportImage(ImageIterator imageUpperLeft, ImageIterator imageLowerRight, ImageAccessor imageAccessor,
     *                          const ImageExportInfo& exportInfo)
     *     }
     * \endcode
     *
     * <B>Usage</B>
     *
     * <B>\#include \<vigra/impex.hxx\></B>
     *
     * Namespace: vigra
     * \code
     *     BRGBImage image(width, height);
     *     ...
     *
     *     // write as JPEG image, using compression quality 80
     *     exportImage(srcImageRange(image),
     *                 ImageExportInfo("my-image.jpg").setCompression("80"));
     *
     *     // Force it to a particular pixel type.  The pixel type must be supported by the
     *     // desired image file format, otherwise an \ref vigra::PreconditionViolation
     *     // exception will be thrown.
     *     exportImage(srcImageRange(image),
     *                 ImageExportInfo("my-INT16-image.tif").setPixelType("INT16"));
     * \endcode
     *
     * <B>Preconditions</B>
     *
     * - The image file must be writable and
     * - the file type must be one of the supported file types.
     */
    doxygen_overloaded_function(template <...> inline void exportImage)


    template <class ImageIterator, class ImageAccessor>
    inline void
    exportImage(ImageIterator image_upper_left, ImageIterator image_lower_right, ImageAccessor image_accessor,
                const ImageExportInfo& export_info)
    {
        typedef typename ImageAccessor::value_type ImageValueType;
        typedef typename NumericTraits<ImageValueType>::isScalar is_scalar;

        try
        {
            detail::exportImage(image_upper_left, image_lower_right, image_accessor,
                        export_info,
                        is_scalar());
        }
        catch (Encoder::TIFFCompressionException&)
        {
            ImageExportInfo info(export_info);

            info.setCompression("");
            detail::exportImage(image_upper_left, image_lower_right, image_accessor,
                                   info,
                                   is_scalar());
        }
    }


    template <class ImageIterator, class ImageAccessor>
    inline void
    exportImage(const vigra::triple<ImageIterator, ImageIterator, ImageAccessor>& image,
                const ImageExportInfo& export_info)
    {
        exportImage(image.first, image.second, image.third,
                    export_info);
    }

/** @} */

} // end namespace vigra

#endif // VIGRA_IMPEX_HXX
