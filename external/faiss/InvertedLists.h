/**
 * Copyright (c) 2015-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD+Patents license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#ifndef FAISS_INVERTEDLISTS_IVF_H
#define FAISS_INVERTEDLISTS_IVF_H

/**
 * Definition of inverted lists + a few common classes that implement
 * the interface.
 */

#include <vector>
#include "Index.h"


namespace faiss {

/** Table of inverted lists
 * multithreading rules:
 * - concurrent read accesses are allowed
 * - concurrent update accesses are allowed
 * - for resize and add_entries, only concurrent access to different lists
 *   are allowed
 */
struct InvertedLists {
    typedef Index::idx_t idx_t;

    size_t nlist;             ///< number of possible key values
    size_t code_size;         ///< code size per vector in bytes

    InvertedLists (size_t nlist, size_t code_size);

    /*************************
     *  Read only functions */

    /// get the size of a list
    virtual size_t list_size(size_t list_no) const = 0;

    /** get the codes for an inverted list
     * must be released by release_codes
     *
     * @return codes    size list_size * code_size
     */
    virtual const uint8_t * get_codes (size_t list_no) const = 0;

    /** get the ids for an inverted list
     * must be released by release_ids
     *
     * @return ids      size list_size
     */
    virtual const idx_t * get_ids (size_t list_no) const = 0;

    /// release codes returned by get_codes (default implementation is nop
    virtual void release_codes (const uint8_t *codes) const;

    /// release ids returned by get_ids
    virtual void release_ids (const idx_t *ids) const;

    /// @return a single id in an inverted list
    virtual idx_t get_single_id (size_t list_no, size_t offset) const;

    /// @return a single code in an inverted list
    /// (should be deallocated with release_codes)
    virtual const uint8_t * get_single_code (
                size_t list_no, size_t offset) const;

    /// prepare the following lists (default does nothing)
    /// a list can be -1 hence the signed long
    virtual void prefetch_lists (const long *list_nos, int nlist) const;

    /*************************
     * writing functions     */

    /// add one entry to an inverted list
    virtual size_t add_entry (size_t list_no, idx_t theid,
                              const uint8_t *code);

    virtual size_t add_entries (
           size_t list_no, size_t n_entry,
           const idx_t* ids, const uint8_t *code) = 0;

    virtual void update_entry (size_t list_no, size_t offset,
                               idx_t id, const uint8_t *code);

    virtual void update_entries (size_t list_no, size_t offset, size_t n_entry,
                                 const idx_t *ids, const uint8_t *code) = 0;

    virtual void resize (size_t list_no, size_t new_size) = 0;

    virtual void reset ();

    /// move all entries from oivf (empty on output)
    void merge_from (InvertedLists *oivf, size_t add_id);

    virtual ~InvertedLists ();

    /**************************************
     * Scoped inverted lists (for automatic deallocation)
     *
     * instead of writing:
     *
     *     uint8_t * codes = invlists->get_codes (10);
     *     ... use codes
     *     invlists->release_codes(codes)
     *
     * write:
     *
     *    ScopedCodes codes (invlists, 10);
     *    ... use codes.get()
     *    // release called automatically when codes goes out of scope
     *
     * the following function call also works:
     *
     *    foo (123, ScopedCodes (invlists, 10).get(), 456);
     *
     */

    struct ScopedIds {
        const InvertedLists *il;
        const idx_t *ids;

        ScopedIds (const InvertedLists *il, size_t list_no):
        il (il), ids (il->get_ids (list_no))
        {}

        const idx_t *get() {return ids; }

        idx_t operator [] (size_t i) const {
            return ids[i];
        }

        ~ScopedIds () {
            il->release_ids (ids);
        }
    };

    struct ScopedCodes {
        const InvertedLists *il;
        const uint8_t *codes;

        ScopedCodes (const InvertedLists *il, size_t list_no):
        il (il), codes (il->get_codes (list_no))
        {}

        ScopedCodes (const InvertedLists *il, size_t list_no, size_t offset):
        il (il), codes (il->get_single_code (list_no, offset))
        {}

        const uint8_t *get() {return codes; }

        ~ScopedCodes () {
            il->release_codes (codes);
        }
    };


};


/// simple (default) implementation as an array of inverted lists
struct ArrayInvertedLists: InvertedLists {
    std::vector < std::vector<uint8_t> > codes; // binary codes, size nlist
    std::vector < std::vector<idx_t> > ids;  ///< Inverted lists for indexes

    ArrayInvertedLists (size_t nlist, size_t code_size);

    size_t list_size(size_t list_no) const override;
    const uint8_t * get_codes (size_t list_no) const override;
    const idx_t * get_ids (size_t list_no) const override;

    size_t add_entries (
           size_t list_no, size_t n_entry,
           const idx_t* ids, const uint8_t *code) override;

    void update_entries (size_t list_no, size_t offset, size_t n_entry,
                         const idx_t *ids, const uint8_t *code) override;

    void resize (size_t list_no, size_t new_size) override;

    virtual ~ArrayInvertedLists ();
};


/// inverted lists built as the concatenation of a set of invlists
/// (read-only)
struct ConcatenatedInvertedLists: InvertedLists {

    std::vector<const InvertedLists *>ils;

    /// build InvertedLists by concatenating nil of them
    ConcatenatedInvertedLists (int nil, const InvertedLists **ils);

    size_t list_size(size_t list_no) const override;
    const uint8_t * get_codes (size_t list_no) const override;
    const idx_t * get_ids (size_t list_no) const override;

    void release_codes (const uint8_t *codes) const override;
    void release_ids (const idx_t *ids) const override;

    idx_t get_single_id (size_t list_no, size_t offset) const override;

    const uint8_t * get_single_code (
           size_t list_no, size_t offset) const override;

    size_t add_entries (
           size_t list_no, size_t n_entry,
           const idx_t* ids, const uint8_t *code) override;

    void update_entries (size_t list_no, size_t offset, size_t n_entry,
                         const idx_t *ids, const uint8_t *code) override;

    void resize (size_t list_no, size_t new_size) override;

};

} // namespace faiss


#endif
