/**
 * Author: Fernando Vilas 2025
 *
 * Borrowed heavily from the nixiz blog implementation and updated cor C++23
 */

#pragma once

#include <condition_variable>
#include <deque>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

/**
 * This interface is similar to the boost::asio interface
 * but using C++23 and not as many features.
 * Whenever the C++ standard finally adds thread_pool then
 * this class should be deleted.
 */

class thread_pool
{
public:
    explicit thread_pool(int num_threads = 1)
    {
        for (int i = 0; i < num_threads; i++) {
            mPool.emplace_back(&thread_pool::run, this);
        }
    }
    ~thread_pool() {
        for (auto& th : mPool) {
            th.request_stop();
        }
        for (auto& th : mPool) {
            th.join();
        }
    }

    thread_pool(const thread_pool&) = delete;
    auto operator=(const thread_pool&) -> thread_pool& = delete;

    thread_pool(thread_pool&&) = delete;
    auto operator=(thread_pool&&) -> thread_pool& = delete;

    void post(std::packaged_task<void()> job) {
        const std::unique_lock lock(mGuard);
        mPendingJobs.emplace_back(std::move(job));
        mCV.notify_one();
    }


    struct use_future_tag {};

    template <class Fn>
    static constexpr auto use_future(Fn&& func) {
        return std::make_tuple(use_future_tag{}, std::forward<Fn>(func));
    }

    template <class Executor, class Fn>
    static void post(Executor& exec, Fn&& func)
    {
        using return_type = decltype(func());
        static_assert(std::is_void_v<return_type>, "posting functions with return types must be used with \"use_future\" tag.");
        std::packaged_task<void()> task(std::forward<Fn>(func));
        exec.post(std::move(task));
    }

    template <class Executor, class Fn>
    [[nodiscard]] static auto
    post(Executor& exec, std::tuple<use_future_tag, Fn>&& tpl) // -> decltype(auto)
    {
        using return_type = std::invoke_result_t<Fn>;
        auto&& [_, func] = tpl;
        if constexpr (std::is_void_v<return_type>)
        {
            std::packaged_task<void()> tsk(std::move(func));
            auto ret_future = tsk.get_future();
            exec.post(std::move(tsk));
            return ret_future;
        }
        else
        {
            struct forwarder_t {
                explicit forwarder_t(Fn&& fn) : tsk(std::forward<Fn>(fn)) {}
                void operator()(std::shared_ptr<std::promise<return_type>> promise) noexcept
                {
                    promise->set_value(tsk());
                }
            private:
                std::decay_t<Fn> tsk;
            } forwarder(std::forward<Fn>(func));

            auto promise = std::make_shared<std::promise<return_type>>();
            auto ret_future = promise->get_future();
            std::packaged_task<void()> tsk([promise = std::move(promise), forwarder = std::move(forwarder)] () mutable {
                forwarder(promise);
            });
            exec.post(std::move(tsk));
            return ret_future;
        }
    }

private:
    void run(std::stop_token stoken) noexcept // NOLINT(performance-unnecessary-value-param)
    {
        while (!stoken.stop_requested())
        {
            thread_local std::packaged_task<void()> job;
            {
                std::unique_lock lock(mGuard);
                mCV.wait(lock, stoken, [&]{ return !mPendingJobs.empty(); });
                if (stoken.stop_requested()) {
                    break;
                }
                job.swap(mPendingJobs.front());
                mPendingJobs.pop_front();
            }
            job();
        }
    }

    std::vector<std::jthread> mPool;
    std::condition_variable_any mCV;
    std::mutex mGuard;
    std::deque<std::packaged_task<void()>> mPendingJobs;
};
